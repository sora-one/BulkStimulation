# ============================================================================
# RNA-Seq数据的混杂因素校正
# ============================================================================
# 本脚本包含用于RNA-seq数据中已知和未知混杂因素的各种校正方法的函数，
# 参照Cote等人(2022, DOI: 10.1186/s13059-022-02606-0)的研究
# ============================================================================

# 加载所需的库
suppressPackageStartupMessages({
  library(DESeq2)       # 用于方差稳定转换
  library(limma)        # 用于removeBatchEffect和voom函数
  library(sva)          # 用于SVA和num.sv函数
  library(matrixStats)  # 用于rowVars函数
  library(PEER)         # 用于PEER校正
  library(RUVcorr)      # 用于RUVCorr校正
  library(BiocParallel) # 用于并行处理
  library(ggplot2)      # 用于可视化
  library(dplyr)        # 用于数据操作
})

#' 对RNA-Seq数据应用多种混杂因素校正方法
#'
#' 此函数对RNA-seq数据应用各种校正方法，以消除技术和生物混杂因素的影响。
#'
#' @param sim_data 包含模拟RNA-seq数据的列表，包含以下元素:
#'   - bulk_counts: 原始基因表达计数矩阵（基因 x 样本）
#'   - bulk_metadata: 样本元数据数据框
#'   - gene_metadata: 基因元数据数据框
#'
#' @return 一个校正后表达矩阵的列表，每个矩阵以应用的校正方法命名。
#'
apply_all_corrections <- function(sim_data) {
  
  # 从输入列表中提取数据
  bulk_counts <- sim_data$bulk_counts
  bulk_metadata <- sim_data$bulk_metadata
  gene_metadata <- sim_data$gene_metadata
  
  # 检查输入数据
  if (is.null(bulk_counts) || is.null(bulk_metadata) || is.null(gene_metadata)) {
    stop("输入数据必须包含bulk_counts、bulk_metadata和gene_metadata")
  }
  
  # 确保样本名称在计数矩阵和元数据之间匹配
  if (!identical(colnames(bulk_counts), rownames(bulk_metadata))) {
    stop("bulk_counts中的样本名称与bulk_metadata中的不匹配")
  }
  
  # 确保基因名称在计数矩阵和元数据之间匹配
  if (!identical(rownames(bulk_counts), rownames(gene_metadata))) {
    stop("bulk_counts中的基因名称与gene_metadata中的不匹配")
  }
  
  # 创建一个列表存储所有校正后的矩阵
  corrected_matrices <- list()
  
  # ============================================================================
  # 1. 数据标准化
  # ============================================================================
  message("执行方差稳定转换...")
  
  # 创建DESeqDataSet对象
  dds <- DESeqDataSetFromMatrix(
    countData = bulk_counts,
    colData = bulk_metadata,
    design = ~ group
  )
  
  # 应用方差稳定转换
  vsd <- vst(dds, blind = FALSE)
  
  # 提取标准化后的矩阵
  normalized_matrix <- assay(vsd)
  
  # 添加到结果列表
  corrected_matrices$normalized_uncorrected <- normalized_matrix
  
  # ============================================================================
  # 2. 已知协变量回归
  # ============================================================================
  message("应用已知协变量回归...")
  
  # 为生物学因素(group)创建设计矩阵
  design <- model.matrix(~ group, data = bulk_metadata)
  
  # 识别细胞类型比例列（以"prop_"开头）
  prop_cols <- grep("^prop_", names(bulk_metadata), value = TRUE)
  
  if (length(prop_cols) == 0) {
    warning("在bulk_metadata中未找到细胞类型比例列")
    covariates <- model.matrix(~ batch, data = bulk_metadata)[, -1, drop = FALSE]
  } else {
    # 为批次和细胞类型比例创建协变量矩阵
    covariates <- model.matrix(
      as.formula(paste("~ batch +", paste(prop_cols, collapse = " + "))),
      data = bulk_metadata
    )[, -1, drop = FALSE]  # 移除截距
  }
  
  # 应用removeBatchEffect校正已知协变量
  known_covariates_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = covariates,
    design = design
  )
  
  # 添加到结果列表
  corrected_matrices$known_covariates <- known_covariates_matrix
  
  # ============================================================================
  # 3. 主成分回归 (PCR)
  # ============================================================================
  message("应用主成分回归...")
  
  # 过滤低方差基因
  gene_vars <- rowVars(normalized_matrix)
  high_var_genes <- which(gene_vars > quantile(gene_vars, 0.5))
  
  # 对高变异基因进行PCA
  pca_data <- prcomp(t(normalized_matrix[high_var_genes, ]), scale = TRUE, center = TRUE)
  
  # 确定要移除的显著PC数量
  # 仅为生物学因素(group)创建设计矩阵
  mod <- model.matrix(~ group, data = bulk_metadata)
  mod0 <- model.matrix(~ 1, data = bulk_metadata)
  
  # 估计替代变量的数量
  n_sv <- num.sv(normalized_matrix, mod, method = "be", B = 20)
  message(paste("估计的显著PC数量:", n_sv))
  
  # 创建要回归的PC矩阵
  pcs_to_regress <- pca_data$x[, 1:n_sv, drop = FALSE]
  
  # 应用removeBatchEffect移除显著PC的影响
  pc_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = pcs_to_regress,
    design = mod
  )
  
  # 添加到结果列表
  corrected_matrices$pc_correction <- pc_corrected_matrix
  
  # ============================================================================
  # 4. 替代变量分析 (SVA)
  # ============================================================================
  message("应用替代变量分析...")
  
  # 执行SVA
  sva_result <- sva(normalized_matrix, mod, mod0, n.sv = n_sv)
  
  # 提取替代变量
  sv <- sva_result$sv
  
  # 应用removeBatchEffect移除替代变量的影响
  sva_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = sv,
    design = mod
  )
  
  # 添加到结果列表
  corrected_matrices$sva_correction <- sva_corrected_matrix
  
  # ============================================================================
  # 5. PEER (表达残差的概率估计)
  # ============================================================================
  message("应用PEER校正...")
  
  # 根据样本大小确定PEER因子的数量
  n_samples <- ncol(bulk_counts)
  if (n_samples < 150) {
    n_factors <- 15
  } else if (n_samples < 250) {
    n_factors <- 30
  } else if (n_samples < 350) {
    n_factors <- 45
  } else {
    n_factors <- 60
  }
  message(paste("对", n_samples, "个样本使用", n_factors, "个PEER因子"))
  
  # 设置PEER模型
  peer_model <- PEER()
  PEER_setPhenoMean(peer_model, t(normalized_matrix))
  
  # 添加已知协变量(group)以保护生物学信号
  PEER_setNk(peer_model, n_factors)
  PEER_update(peer_model)
  
  # 获取残差（校正后的数据）
  residuals <- t(PEER_getResiduals(peer_model))
  rownames(residuals) <- rownames(normalized_matrix)
  colnames(residuals) <- colnames(normalized_matrix)
  
  # 添加到结果列表
  corrected_matrices$peer_correction <- residuals
  
  # ============================================================================
  # 6. RUVCorr
  # ============================================================================
  message("应用RUVCorr校正...")
  
  # 寻找阴性对照基因（不在任何模块中的基因）
  module_cols <- grep("^in_module_", colnames(gene_metadata), value = TRUE)
  
  if (length(module_cols) == 0) {
    warning("在gene_metadata中未找到in_module列；对RUVCorr使用所有基因")
    control_genes <- rep(TRUE, nrow(gene_metadata))
  } else {
    # 选择不在任何模块中的基因
    control_genes <- rowSums(gene_metadata[, module_cols, drop = FALSE]) == 0
  }
  
  # 如果未找到阴性对照基因，则使用方差最低的基因
  if (sum(control_genes) < 10) {
    warning("找到的阴性对照基因少于10个；使用方差最低的10%的基因")
    control_genes <- gene_vars <= quantile(gene_vars, 0.1)
  }
  
  message(paste("为RUVCorr使用", sum(control_genes), "个阴性对照基因"))
  
  # 设置要移除的非期望变异因子数量(k)
  k <- min(5, ncol(bulk_counts) / 10)  # 使用5或样本大小的10%，取较小值
  
  # 应用RUVCorr
  ruv_result <- RUVcorr(Y = normalized_matrix, 
                        X = model.matrix(~ group, data = bulk_metadata),
                        ctl = control_genes, 
                        k = k)
  
  # 提取校正后的矩阵
  ruv_corrected_matrix <- ruv_result$Y_corr
  
  # 添加到结果列表
  corrected_matrices$ruv_corr <- ruv_corrected_matrix
  
  # ============================================================================
  # 返回所有校正后的矩阵
  # ============================================================================
  message("所有校正方法已成功应用。")
  return(corrected_matrices)
}

# ============================================================================
# 使用示例
# ============================================================================

#' apply_all_corrections函数的使用示例
example_usage <- function() {
  # 加载模拟数据
  # 替换为你的实际文件路径
  sim_data <- readRDS("path/to/your/simulated_data.rds")
  
  # 应用所有校正方法
  corrected_data <- apply_all_corrections(sim_data)
  
  # 打印每个校正矩阵的维度
  for (method_name in names(corrected_data)) {
    matrix_dim <- dim(corrected_data[[method_name]])
    cat(sprintf("%s: %d基因 × %d样本\n", 
                method_name, matrix_dim[1], matrix_dim[2]))
  }
  
  # 返回校正后的数据
  return(corrected_data)
}

# 取消注释以运行示例
# results <- example_usage()
