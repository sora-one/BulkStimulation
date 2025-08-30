# ============================================================================
# RNA-Seq数据的混杂因素校正
# ============================================================================
# 本脚本包含用于RNA-seq数据中已知和未知混杂因素的各种校正方法的函数，
# 参照Cote等人(2022, DOI: 10.1186/s13059-022-02606-0)的研究
# ============================================================================
# 主要功能：apply_all_corrections
# 作用：对输入的RNA-seq表达数据应用多种混杂因素校正方法，返回每种方法校正后的表达矩阵
# 方法包括：方差稳定转换、已知协变量回归、主成分回归(PCR)、替代变量分析(SVA)、PEER、RUVCorr
# ============================================================================
# 所需R包：DESeq2, limma, sva, matrixStats, PEER, RUVcorr, BiocParallel, ggplot2, dplyr
# ============================================================================

# 加载所需的库，避免启动消息
suppressPackageStartupMessages({
  library(DESeq2)       # 用于RNA-seq计数数据的方差稳定转换 (vst)
  library(limma)        # 用于批次效应移除 (removeBatchEffect) 和线性建模
  library(sva)          # 用于替代变量分析 (SVA) 及估算显著成分 (num.sv)
  library(matrixStats)  # 用于高效计算基因表达矩阵的行方差 (rowVars)
  library(PEER)         # 用于PEER模型校正表达数据中的隐含混杂因素
  library(RUVcorr)      # 用于RUVCorr方法，利用阴性对照基因校正
  library(BiocParallel) # 用于并行计算环境
  library(ggplot2)      # 用于可视化分析（未在主流程中使用，准备可选绘图）
  library(dplyr)        # 用于数据处理和变换
})

#' 对RNA-Seq数据应用多种混杂因素校正方法
#'
#' @param sim_data 包含模拟RNA-seq数据的列表，需包含以下元素:
#'   - bulk_counts: 原始基因表达计数矩阵（基因 x 样本，行为基因，列为样本）
#'   - bulk_metadata: 样本元数据（行名为样本，需含"group"和"batch"等信息）
#'   - gene_metadata: 基因元数据（行名为基因，需含模块信息等）
#'
#' @return 一个校正后表达矩阵的列表，每个矩阵以所应用的校正方法命名。
#'         包括未校正标准化(normalized_uncorrected)、已知协变量(known_covariates)、主成分回归(pc_correction)、SVA(sva_correction)、PEER(peer_correction)、RUVCorr(ruv_corr)
#'
apply_all_corrections <- function(sim_data) {
  
  # ============================== 输入数据检查 ================================
  # 提取计数矩阵与元数据
  bulk_counts <- sim_data$bulk_counts
  bulk_metadata <- sim_data$bulk_metadata
  gene_metadata <- sim_data$gene_metadata
  
  # 检查输入完整性
  if (is.null(bulk_counts) || is.null(bulk_metadata) || is.null(gene_metadata)) {
    stop("输入数据必须包含bulk_counts、bulk_metadata和gene_metadata")
  }
  
  # 检查样本名称是否一致
  if (!identical(colnames(bulk_counts), rownames(bulk_metadata))) {
    stop("bulk_counts中的样本名称与bulk_metadata中的不匹配")
  }
  
  # 检查基因名称是否一致
  if (!identical(rownames(bulk_counts), rownames(gene_metadata))) {
    stop("bulk_counts中的基因名称与gene_metadata中的不匹配")
  }
  
  # 初始化用于存储所有校正后表达矩阵的列表
  corrected_matrices <- list()
  
  # ============================================================================
  # 1. 数据标准化 (方差稳定转换)
  # ============================================================================
  # 使用DESeq2的vst对原始计数数据进行标准化，消除测序深度等技术噪音
  message("执行方差稳定转换...")
  
  # 构造DESeqDataSet对象。design = ~ group 表示以"group"为生物学变量
  dds <- DESeqDataSetFromMatrix(
    countData = bulk_counts,
    colData = bulk_metadata,
    design = ~ group
  )
  
  # 应用方差稳定转换 (vst)，blind=FALSE表示考虑design中的变量
  vsd <- vst(dds, blind = FALSE)
  
  # 提取标准化后的表达矩阵
  normalized_matrix <- assay(vsd)
  
  # 保存未校正的标准化表达矩阵
  corrected_matrices$normalized_uncorrected <- normalized_matrix
  
  # ============================================================================
  # 2. 已知协变量回归 (批次/细胞类型比例)
  # ============================================================================
  # 利用limma::removeBatchEffect校正已知技术/生物学协变量
  message("应用已知协变量回归...")
  
  # 构建设计矩阵，包含生物学因素 group
  design <- model.matrix(~ group, data = bulk_metadata)
  
  # 搜索样本元数据中所有“prop_”开头的细胞类型比例列
  prop_cols <- grep("^prop_", names(bulk_metadata), value = TRUE)
  
  # 构建协变量矩阵：如有prop_列，则包含batch与prop_，否则仅batch
  if (length(prop_cols) == 0) {
    warning("在bulk_metadata中未找到细胞类型比例列")
    covariates <- model.matrix(~ batch, data = bulk_metadata)[, -1, drop = FALSE]
  } else {
    covariates <- model.matrix(
      as.formula(paste("~ batch +", paste(prop_cols, collapse = " + "))),
      data = bulk_metadata
    )[, -1, drop = FALSE]  # 移除截距
  }
  
  # 用limma校正已知协变量影响
  known_covariates_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = covariates,
    design = design
  )
  
  # 保存校正结果
  corrected_matrices$known_covariates <- known_covariates_matrix
  
  # ============================================================================
  # 3. 主成分回归 (Principal Component Regression, PCR)
  # ============================================================================
  # 用主成分回归方法校正高变异的隐含混杂因素
  message("应用主成分回归...")
  
  # 计算基因表达矩阵的行方差，筛选高变异基因
  gene_vars <- rowVars(normalized_matrix)
  high_var_genes <- which(gene_vars > quantile(gene_vars, 0.5))
  
  # 对高变异基因的表达数据做PCA，提取主成分
  pca_data <- prcomp(t(normalized_matrix[high_var_genes, ]), scale = TRUE, center = TRUE)
  
  # 构建设计矩阵 mod（含group），mod0（仅截距）
  mod <- model.matrix(~ group, data = bulk_metadata)
  mod0 <- model.matrix(~ 1, data = bulk_metadata)
  
  # 用sva::num.sv估算需要移除的显著主成分数量
  n_sv <- num.sv(normalized_matrix, mod, method = "be", B = 20)
  message(paste("估计的显著PC数量:", n_sv))
  
  # 取前n_sv个主成分，作为待回归的隐含混杂因素
  pcs_to_regress <- pca_data$x[, 1:n_sv, drop = FALSE]
  
  # 用limma校正主成分影响
  pc_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = pcs_to_regress,
    design = mod
  )
  
  # 保存PCR校正结果
  corrected_matrices$pc_correction <- pc_corrected_matrix
  
  # ============================================================================
  # 4. 替代变量分析 (Surrogate Variable Analysis, SVA)
  # ============================================================================
  # 自动识别并校正表达数据中的隐含混杂因素
  message("应用替代变量分析...")
  
  # 执行SVA分析，提取替代变量
  sva_result <- sva(normalized_matrix, mod, mod0, n.sv = n_sv)
  sv <- sva_result$sv
  
  # 用limma校正替代变量影响
  sva_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = sv,
    design = mod
  )
  
  # 保存SVA校正结果
  corrected_matrices$sva_correction <- sva_corrected_matrix
  
  # ============================================================================
  # 5. PEER (Probabilistic Estimation of Expression Residuals)
  # ============================================================================
  # 用PEER模型校正表达数据中的隐含混杂因素
  message("应用PEER校正...")
  
  # 根据样本数自动确定PEER因子数量，防止过拟合或欠拟合
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
  
  # 构建PEER模型对象，输入标准化表达矩阵
  peer_model <- PEER()
  PEER_setPhenoMean(peer_model, t(normalized_matrix))
  
  # 设置PEER因子数量
  PEER_setNk(peer_model, n_factors)
  # PEER可接受已知协变量（如group），但此处未显式设置
  
  # 运行PEER模型，估算隐含因子
  PEER_update(peer_model)
  
  # 提取PEER校正后的残差矩阵，作为去除隐含混杂因素后的表达数据
  residuals <- t(PEER_getResiduals(peer_model))
  rownames(residuals) <- rownames(normalized_matrix)
  colnames(residuals) <- colnames(normalized_matrix)
  
  # 保存PEER校正结果
  corrected_matrices$peer_correction <- residuals
  
  # ============================================================================
  # 6. RUVCorr (Removal of Unwanted Variation using Negative Control Genes)
  # ============================================================================
  # 用阴性对照基因（未参与生物学模块）校正非期望变异
  message("应用RUVCorr校正...")
  
  # 搜索基因元数据中的模块注释列，以"in_module_"开头
  module_cols <- grep("^in_module_", colnames(gene_metadata), value = TRUE)
  
  # 若未找到模块注释，则所有基因作为阴性对照
  if (length(module_cols) == 0) {
    warning("在gene_metadata中未找到in_module列；对RUVCorr使用所有基因")
    control_genes <- rep(TRUE, nrow(gene_metadata))
  } else {
    # 选择未在任何模块中的基因作为阴性对照（各模块标志列都为0）
    control_genes <- rowSums(gene_metadata[, module_cols, drop = FALSE]) == 0
  }
  
  # 若阴性对照基因数量太少，选方差最低的10%基因作为替代
  if (sum(control_genes) < 10) {
    warning("找到的阴性对照基因少于10个；使用方差最低的10%的基因")
    control_genes <- gene_vars <= quantile(gene_vars, 0.1)
  }
  
  message(paste("为RUVCorr使用", sum(control_genes), "个阴性对照基因"))
  
  # 设置RUVCorr要移除的非期望变异因子数量k（最多5或样本数的10%）
  k <- min(5, ncol(bulk_counts) / 10)
  
  # 应用RUVCorr校正表达矩阵
  ruv_result <- RUVcorr(Y = normalized_matrix, 
                        X = model.matrix(~ group, data = bulk_metadata),
                        ctl = control_genes, 
                        k = k)
  
  # 提取校正后的表达矩阵
  ruv_corrected_matrix <- ruv_result$Y_corr
  
  # 保存RUVCorr校正结果
  corrected_matrices$ruv_corr <- ruv_corrected_matrix
  
  # ============================================================================
  # 返回所有校正后表达矩阵
  # ============================================================================
  message("所有校正方法已成功应用。")
  return(corrected_matrices)
}

# ============================================================================
# 使用示例
# ============================================================================
# 此函数演示如何加载模拟数据，应用所有校正方法，并输出每种方法校正后矩阵的维度
# 实际使用时请替换数据路径
# ============================================================================
#' apply_all_corrections函数的使用示例
example_usage <- function() {
  # 加载模拟数据
  # 请将"path/to/your/simulated_data.rds"替换为实际数据文件路径
  sim_data <- readRDS("path/to/your/simulated_data.rds")
  
  # 对模拟数据应用所有混杂因素校正方法
  corrected_data <- apply_all_corrections(sim_data)
  
  # 遍历所有校正结果，打印每个校正矩阵的维度信息
  for (method_name in names(corrected_data)) {
    matrix_dim <- dim(corrected_data[[method_name]])
    cat(sprintf("%s: %d基因 × %d样本\n", 
                method_name, matrix_dim[1], matrix_dim[2]))
  }
  
  # 返回校正后的所有表达矩阵
  return(corrected_data)
}

# 如需运行示例，请取消下行注释
# results <- example_usage()
