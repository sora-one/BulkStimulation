#' @title RNA-seq 数据模拟管道
#' @description 一个用于生成模拟单细胞和bulk RNA-seq数据的综合管道
#' @author XushengWan
#' @date 2025-08-25
#' 
#' 本脚本提供了一个完整的框架，用于:
#' 1. 生成具有预定义细胞类型和基因共表达模块的单细胞数据
#' 2. 根据不同的混合方案将单细胞数据聚合为bulk RNA-seq样本
#' 3. 支持系统性参数探索和并行计算
#'

# 加载必要的库 -------------------------------------------------------------
suppressPackageStartupMessages({
  library(BiocParallel)       # 用于并行计算
  library(splatter)           # 用于单细胞RNA-seq模拟
  library(scDesign3)          # 用于增强共表达网络
  library(SingleCellExperiment) # 用于单细胞数据处理
  library(scater)             # 用于单细胞数据可视化和预处理
  library(Matrix)             # 用于高效处理稀疏矩阵
  library(ggplot2)            # 用于数据可视化
  library(dplyr)              # 用于数据操作
  library(DESeq2)             # 用于后续差异表达分析
})

# 辅助函数: 日志打印 -------------------------------------------------------
#' @description 打印带有时间戳的日志信息
log_message <- function(msg, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  # 处理格式化字符串
  if (...length() > 0) {
    msg <- sprintf(msg, ...)
  }
  cat(sprintf("[%s] [INFO] %s\n", timestamp, msg))
}

# 模块1: 参数配置模块 ------------------------------------------------------
#' @description 配置单细胞RNA-seq模拟的参数
#' @param intra_group_distance 组内细胞距离 (0.05-0.5)
#' @param inter_group_distance 组间细胞距离 (0.5-5.0)
#' @param overall_cell_spread 整体细胞分布散度 (0.1-0.5)
#' @param group_separation_strength 组间分离强度 (0.5-2.0)
#' @param de_gene_proportion 差异表达基因比例 (0.1-0.4)
#' @param n_cells 细胞总数
#' @param n_genes 基因总数
#' @param n_groups 细胞组数
#' @param group_probabilities 各组细胞比例
#' @param seed 随机种子
#' @param verbose 是否打印详细信息
#' @return 配置好的参数列表
configure_simulation_params <- function(
    # 细胞间距离控制参数
  intra_group_distance = 0.1,      # 组内细胞间平均距离 (0.05-0.5)
  inter_group_distance = 1.0,      # 组间细胞平均距离 (0.5-5.0)
  overall_cell_spread = 0.2,       # 整体细胞分布散度 (0.1-0.5)
  
  # 差异表达控制参数
  group_separation_strength = 1.0, # 组间分离强度 (0.5-2.0)
  de_gene_proportion = 0.2,        # 差异表达基因比例 (0.1-0.4)
  
  # 基本模拟参数
  n_cells = 3000,                  # 细胞总数
  n_genes = 3000,                  # 基因总数
  n_groups = 6,                    # 细胞组数
  group_probabilities = NULL,      # 各组细胞比例 (如不提供则自动生成)
  
  # 技术噪音控制
  biological_cv = NULL,            # 生物学变异系数 (自动根据距离参数计算)
  library_size_variation = NULL,   # 文库大小变异 (自动根据距离参数计算)
  
  # 随机种子
  seed = 42,
  
  # 输出详细信息
  verbose = TRUE
) {
  # 设置随机种子以确保可重复性
  set.seed(seed)
  
  # 如果没有提供组概率，则生成均衡的概率
  if (is.null(group_probabilities)) {
    group_probabilities <- rep(1/n_groups, n_groups)
  }
  
  # 参数验证
  if (length(group_probabilities) != n_groups) {
    stop("group_probabilities 长度必须等于 n_groups")
  }
  
  if (abs(sum(group_probabilities) - 1) > 0.01) {
    stop("group_probabilities 总和必须等于 1")
  }
  
  # 根据距离参数自动计算技术参数
  if (is.null(biological_cv)) {
    # 组内距离越小，生物学变异越小
    biological_cv <- 0.05 + (intra_group_distance * 0.5)
  }
  
  if (is.null(library_size_variation)) {
    # 整体散度越小，文库大小变异越小
    library_size_variation <- 0.1 + (overall_cell_spread * 0.4)
  }
  
  # 计算差异表达相关参数
  # 组间距离越大，差异表达强度越大
  de_facLoc <- 0.5 + (inter_group_distance * group_separation_strength * 0.4)
  
  # 组内距离影响组内差异表达概率的变异
  de_prob_base <- de_gene_proportion
  de_prob_variation <- intra_group_distance * 0.1
  de_probs <- pmax(0.05, pmin(0.4,
                              de_prob_base + rnorm(n_groups, 0, de_prob_variation)))
  
  # 创建 Splatter 参数对象
  params <- newSplatParams(
    batchCells = n_cells,
    nGenes = n_genes,
    group.prob = group_probabilities,
    de.prob = de_probs,
    de.facLoc = de_facLoc,
    seed = seed,
    bcv.common = biological_cv,
    bcv.df = 60,
    lib.loc = 11,
    lib.scale = library_size_variation
  )
  
  # 输出参数信息
  if (verbose) {
    log_message("=== 单细胞模拟参数配置完成 ===")
    log_message("细胞间距离控制: 组内=%.2f, 组间=%.2f, 整体散度=%.2f", 
                intra_group_distance, inter_group_distance, overall_cell_spread)
    log_message("差异表达参数: 组间分离强度=%.2f, DE基因比例=%.2f, 计算的DE强度=%.3f", 
                group_separation_strength, de_gene_proportion, de_facLoc)
    log_message("技术参数: 生物变异系数=%.3f, 文库大小变异=%.3f", 
                biological_cv, library_size_variation)
    log_message("基本设置: 细胞数=%d, 基因数=%d, 组数=%d", 
                n_cells, n_genes, n_groups)
  }
  
  # 返回配置好的参数对象和额外信息
  result <- list(
    splatter_params = params,
    distance_config = list(
      intra_group_distance = intra_group_distance,
      inter_group_distance = inter_group_distance,
      overall_cell_spread = overall_cell_spread,
      group_separation_strength = group_separation_strength
    ),
    computed_params = list(
      biological_cv = biological_cv,
      library_size_variation = library_size_variation,
      de_facLoc = de_facLoc,
      de_probs = de_probs
    )
  )
  
  return(result)
}

# 模块2: 单细胞数据生成模块 -------------------------------------------------
#' @description 生成模拟的单细胞RNA-seq数据
#' @param params 由configure_simulation_params返回的参数列表
#' @param co_expression_modules 共表达模块的列表
#' @param seed 随机种子
#' @param verbose 是否打印详细信息
#' @return SingleCellExperiment对象
generate_single_cell_data <- function(
    params,
    co_expression_modules = NULL,
    seed = 42,
    verbose = TRUE
) {
  # 设置随机种子
  set.seed(seed)
  
  if (verbose) log_message("开始生成单细胞数据...")
  
  # 使用Splatter生成基本的单细胞数据
  sim <- splatSimulate(
    params$splatter_params,
    method = "groups",
    verbose = FALSE
  )
  
  # 确保细胞类型标签正确设置
  n_groups <- length(params$splatter_params@group.prob)
  cell_type_levels <- paste0("Group", 1:n_groups)
  
  # 设置细胞类型标签
  if (all(grepl("^Group", sim$Group))) {
    colData(sim)$CellType <- factor(sim$Group, levels = cell_type_levels)
  } else if (is.numeric(sim$Group)) {
    colData(sim)$CellType <- factor(sim$Group, 
                                    levels = 1:n_groups, 
                                    labels = cell_type_levels)
  } else {
    stop("未知的Group格式!")
  }
  
  # 记录生成的单细胞数据基本信息
  if (verbose) {
    log_message("已生成%d个细胞, %d个基因的单细胞数据", 
                ncol(sim), nrow(sim))
    log_message("细胞类型分布:")
    print(table(sim$CellType))
    log_message("基因名称示例: %s, %s, %s", 
                rownames(sim)[1], rownames(sim)[2], rownames(sim)[3])
  }
  
  # 添加标准化计数和降维
  sim <- logNormCounts(sim)
  sim <- runPCA(sim)
  
  # 标记共表达模块中的基因
  if (!is.null(co_expression_modules)) {
    # 创建一个数据框来存储基因元数据
    gene_metadata <- data.frame(
      gene_id = rownames(sim),
      stringsAsFactors = FALSE
    )
    
    # 标记每个模块中的基因
    for (i in seq_along(co_expression_modules)) {
      module_name <- names(co_expression_modules)[i]
      if (is.null(module_name) || module_name == "") {
        module_name <- paste0("module", i)
      }
      
      module_genes <- co_expression_modules[[i]]
      
      # 验证模块基因是否存在于数据中
      valid_genes <- intersect(module_genes, rownames(sim))
      
      if (length(valid_genes) == 0) {
        warning("模块 '", module_name, "' 中没有有效基因!")
        # 尝试检查是否是格式问题
        if (verbose) {
          log_message("尝试检查基因名格式不匹配问题...")
        }
      }
      
      # 在基因元数据中标记
      gene_metadata[[paste0("in_module_", module_name)]] <- 
        gene_metadata$gene_id %in% valid_genes
      
      if (verbose) {
        log_message("模块'%s'包含%d个有效基因", 
                    module_name, length(valid_genes))
        if (length(valid_genes) > 0) {
          log_message("示例模块基因: %s, %s", 
                      valid_genes[1], 
                      ifelse(length(valid_genes) > 1, valid_genes[2], ""))
        }
      }
    }
    
    # 将基因元数据添加到SingleCellExperiment对象
    rowData(sim) <- gene_metadata
  }
  
  return(sim)
}

# 模块3: 共表达增强模块 ----------------------------------------------------
#' @description 使用scDesign3增强共表达网络结构
#' @param sce 单细胞实验对象
#' @param celltype_col 细胞类型列名
#' @param assay_use 使用的表达矩阵
#' @param co_expression_strength 共表达强度 (0.1-1.0)
#' @param BPPARAM 并行计算参数
#' @param verbose 是否打印详细信息
#' @return 增强后的单细胞数据
enhance_co_expression <- function(
    sce,
    celltype_col = "CellType",
    assay_use = "counts",
    co_expression_strength = 0.7,
    BPPARAM = SerialParam(),
    verbose = TRUE
) {
  if (verbose) log_message("开始增强共表达网络结构...")
  
  # 验证celltype列是否存在
  if (!(celltype_col %in% colnames(colData(sce)))) {
    stop(sprintf("无法找到细胞类型列 '%s'", celltype_col))
  }
  
  # 验证是否有足够的细胞数量用于拟合
  cell_types <- table(sce[[celltype_col]])
  min_cells <- min(cell_types)
  if (min_cells < 10) {
    warning(sprintf("某些细胞类型的细胞数量很少 (最小: %d)，这可能导致拟合不稳定", min_cells))
  }
  
  # 计算要保留多少高变异性基因来调整共表达强度
  n_genes <- nrow(sce)
  genes_to_keep <- max(50, min(round(n_genes * co_expression_strength), n_genes))
  
  # 预处理：计算基因变异性并进行过滤
  gene_var <- rowVars(logcounts(sce))
  top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:genes_to_keep]
  sce_filtered <- sce[top_var_genes,]
  
  if (verbose) {
    log_message("使用%d个高变异基因进行共表达建模 (%.1f%%)", 
                length(top_var_genes), 100 * length(top_var_genes) / n_genes)
  }
  
  # 运行scDesign3，使用适当的参数以减少警告
  suppressWarnings({
    tryCatch({
      scdesign3_result <- scdesign3(
        sce = sce_filtered,            # 使用过滤后的数据
        assay_use = assay_use,
        celltype = celltype_col,      
        pseudotime = NULL,
        spatial = NULL,
        other_covariates = NULL,
        mu_formula = celltype_col,    
        sigma_formula = celltype_col,
        
        # 关键参数调整，解决"Fitting terminated"警告
        family_use = "nb",            # 负二项分布
        n_cores = 1,                  # 单核心内部拟合，避免并行问题
        usebam = FALSE,
        corr_formula = celltype_col,  # 基于细胞类型的相关性
        copula = "gaussian",          # 高斯copula
        DT = FALSE,
        
        # 增加这些参数以提高模型拟合稳定性
        pseudo_obs = TRUE,            # 使用伪观测值提高稳定性
        return_model = FALSE,         # 不返回模型对象节省内存
        nonzerovar = TRUE,            # 只使用非零方差基因
        nonnegative = TRUE,           # 强制非负值
        
        # 并行计算参数不直接传递给内部函数
        BPPARAM = NULL                # 避免BPPARAM不是slot的错误
      )
    }, error = function(e) {
      # 如果scDesign3失败，返回一个简化的结果
      log_message("scDesign3模型拟合失败: %s\n使用原始数据作为替代", e$message)
      
      # 创建一个类似scDesign3结果的简化结构
      list(
        new_count = counts(sce_filtered),
        new_covariate = data.frame(colData(sce_filtered))
      )
    })
  })
  
  if (verbose) log_message("共表达网络增强完成")
  
  # 提取生成的count矩阵
  sc_counts <- scdesign3_result$new_count
  
  # 将增强后的基因与原始基因集合并
  full_counts <- matrix(0, nrow = nrow(sce), ncol = ncol(sce))
  rownames(full_counts) <- rownames(sce)
  colnames(full_counts) <- colnames(sce)
  
  # 复制增强后的基因表达
  common_genes <- intersect(rownames(sc_counts), rownames(full_counts))
  full_counts[common_genes,] <- as.matrix(sc_counts[common_genes,])
  
  # 对于未被增强的基因，保留原始表达值
  missing_genes <- setdiff(rownames(full_counts), rownames(sc_counts))
  if (length(missing_genes) > 0) {
    full_counts[missing_genes,] <- as.matrix(counts(sce)[missing_genes,])
    if (verbose) {
      log_message("%d个低变异性基因保留了原始表达", length(missing_genes))
    }
  }
  
  # 创建新的SingleCellExperiment对象
  enhanced_sce <- SingleCellExperiment(
    assays = list(counts = full_counts),
    colData = colData(sce)
  )
  
  # 添加logcounts
  logcounts(enhanced_sce) <- log1p(counts(enhanced_sce))
  
  # 保留原始数据的行数据
  if (!is.null(rowData(sce)) && ncol(rowData(sce)) > 0) {
    rowData(enhanced_sce) <- rowData(sce)
  }
  
  # 运行PCA
  enhanced_sce <- runPCA(enhanced_sce)
  
  # 简单检查以确认共表达结构
  if (verbose) {
    log_message("增强后的数据包含%d个细胞和%d个基因", 
                ncol(enhanced_sce), nrow(enhanced_sce))
    log_message("细胞类型分布:")
    print(table(enhanced_sce[[celltype_col]]))
  }
  
  return(enhanced_sce)
}

# 模块4: Bulk样本聚合模块 --------------------------------------------------
#' @description 将单细胞数据聚合为bulk RNA-seq样本
#' @param sce 增强后的单细胞数据
#' @param n_bulk_samples 要生成的bulk样本数量
#' @param confounding_level 混杂水平 ("low", "medium", "high")
#' @param n_batches 批次数量
#' @param batch_effect_strength 批次效应强度 (0.1-1.0)
#' @param celltype_col 细胞类型列名
#' @param seed 随机种子
#' @param verbose 是否打印详细信息
#' @return 包含bulk数据和元数据的列表
aggregate_to_bulk <- function(
    sce,
    n_bulk_samples = 30,
    confounding_level = "medium",
    n_batches = 3,
    batch_effect_strength = 0.6,
    celltype_col = "CellType",
    seed = 42,
    verbose = TRUE
) {
  # 设置随机种子
  set.seed(seed)
  
  if (verbose) log_message("开始聚合为bulk样本...")
  
  # 获取细胞类型信息
  cell_types <- unique(sce[[celltype_col]])
  n_cell_types <- length(cell_types)
  
  # 检查细胞类型数量
  if (n_cell_types < 2) {
    stop("需要至少两种细胞类型才能模拟不同的组成比例")
  }
  
  # 定义样本分组（两个主要分组: Group1、Group2）
  n_groups <- 2
  samples_per_group <- n_bulk_samples / n_groups
  
  # 创建平衡的样本信息表
  sample_info <- data.frame(
    sample_id = paste0("Sample_", 1:n_bulk_samples),
    group = rep(paste0("Group", 1:n_groups), each = samples_per_group),
    batch = rep(paste0("Batch", 1:n_batches), length.out = n_bulk_samples),
    stringsAsFactors = FALSE
  )
  
  # 验证批次和组别的平衡性
  if (verbose) {
    log_message("批次-组别交叉表:")
    print(table(sample_info$batch, sample_info$group))
  }
  
  # 初始化bulk样本矩阵
  bulk_counts <- matrix(0, nrow = nrow(sce), ncol = n_bulk_samples)
  rownames(bulk_counts) <- rownames(sce)
  colnames(bulk_counts) <- sample_info$sample_id
  
  # 初始化存储每个样本中细胞类型比例的矩阵
  cell_proportions <- matrix(0, nrow = n_bulk_samples, ncol = n_cell_types)
  colnames(cell_proportions) <- cell_types
  rownames(cell_proportions) <- sample_info$sample_id
  
  # 每个bulk样本包含的细胞数量
  cells_per_sample <- 100
  
  # 根据混杂水平设置细胞类型混合比例
  # 混杂水平定义:
  # - low: 组间细胞类型比例差异大 (较低混杂)
  # - medium: 组间有适度差异
  # - high: 组间差异小 (高度混杂)
  cell_type_bias <- switch(
    confounding_level,
    "low" = 0.5,     # 低混杂 = 高细胞类型组成差异
    "medium" = 0.3,  # 中混杂 = 中等细胞类型组成差异
    "high" = 0.1,    # 高混杂 = 低细胞类型组成差异
    0.3              # 默认为中等
  )
  
  # 为每个组创建细胞类型比例模板
  group_templates <- list()
  for (g in 1:n_groups) {
    # 生成随机的基本比例
    base_props <- runif(n_cell_types)
    base_props <- base_props / sum(base_props)
    
    # 加强某些细胞类型在特定组中的表达
    emphasized_types <- sample(1:n_cell_types, ceiling(n_cell_types/2))
    
    # 根据混杂水平调整比例
    props <- base_props
    props[emphasized_types] <- props[emphasized_types] * (1 + cell_type_bias)
    props[-emphasized_types] <- props[-emphasized_types] * (1 - cell_type_bias)
    props <- props / sum(props)
    
    group_templates[[g]] <- props
  }
  
  # 为每个样本生成bulk数据
  for (i in 1:n_bulk_samples) {
    if (verbose && i %% 10 == 0) {
      log_message("正在处理第%d个bulk样本...", i)
    }
    
    # 确定样本组并获取相应的模板
    group_idx <- ifelse(sample_info$group[i] == "Group1", 1, 2)
    template <- group_templates[[group_idx]]
    
    # 添加小的随机变异到模板
    props <- template + rnorm(n_cell_types, 0, 0.05)
    props <- pmax(0.05, props)  # 确保至少有5%
    props <- props / sum(props) # 归一化
    
    # 存储真实比例
    cell_proportions[i, ] <- props
    
    # 根据比例选择细胞
    selected_cells <- c()
    for (j in 1:n_cell_types) {
      cell_type <- cell_types[j]
      n_cells_to_select <- round(cells_per_sample * props[j])
      
      # 获取该细胞类型的所有细胞
      cells_of_type <- which(sce[[celltype_col]] == cell_type)
      
      # 如果该类型细胞不足，则用有放回抽样
      if (length(cells_of_type) < n_cells_to_select) {
        selected <- sample(cells_of_type, n_cells_to_select, replace = TRUE)
      } else {
        selected <- sample(cells_of_type, n_cells_to_select, replace = FALSE)
      }
      
      selected_cells <- c(selected_cells, selected)
    }
    
    # 汇总选中细胞的表达量
    bulk_counts[, i] <- rowSums(counts(sce)[, selected_cells])
  }
  
  # 添加批次效应
  # 1. 选择批次敏感基因
  n_genes <- nrow(bulk_counts)
  batch_sensitive_genes <- sample(1:n_genes, round(n_genes * 0.6))
  
  # 2. 为每个批次创建特异性基因表达模式
  genes_per_batch <- split(sample(batch_sensitive_genes), 
                           rep(1:n_batches, length.out = length(batch_sensitive_genes)))
  
  # 3. 应用批次效应
  batch_effect_intensity <- 0.5 + batch_effect_strength * 0.5  # 0.75 to 1.25 range
  
  for (i in 1:ncol(bulk_counts)) {
    batch <- sample_info$batch[i]
    batch_num <- as.integer(sub("Batch", "", batch))
    
    # 对该批次的特异性基因应用效应
    effect_genes <- genes_per_batch[[batch_num]]
    bulk_counts[effect_genes, i] <- bulk_counts[effect_genes, i] * 
      rnorm(length(effect_genes), 
            mean = batch_effect_intensity, 
            sd = 0.05 * batch_effect_strength)
    
    # 对其他批次的基因施加抑制效应
    for (j in setdiff(1:n_batches, batch_num)) {
      other_batch_genes <- genes_per_batch[[j]]
      bulk_counts[other_batch_genes, i] <- bulk_counts[other_batch_genes, i] * 
        rnorm(length(other_batch_genes), 
              mean = 1 - 0.4 * batch_effect_strength, 
              sd = 0.05 * batch_effect_strength)
    }
  }
  
  # 全局批次效应：整体文库大小调整
  global_batch_effects <- sapply(1:n_batches, function(b) {
    rnorm(1, mean = 1, sd = 0.2 * batch_effect_strength)
  })
  
  for (i in 1:ncol(bulk_counts)) {
    batch_num <- as.integer(sub("Batch", "", sample_info$batch[i]))
    effect <- global_batch_effects[batch_num]
    bulk_counts[, i] <- bulk_counts[, i] * effect
  }
  
  # 确保count数据为正整数
  bulk_counts[bulk_counts < 0] <- 0
  bulk_counts <- round(bulk_counts)
  
  # 添加质量控制统计
  sample_info$library_size <- colSums(bulk_counts)
  sample_info$detected_genes <- colSums(bulk_counts > 0)
  
  # 添加真实的细胞类型比例到样本信息
  for (j in 1:n_cell_types) {
    sample_info[paste0("prop_", cell_types[j])] <- cell_proportions[, j]
  }
  
  if (verbose) {
    log_message("已生成%d个bulk样本, %d个基因", 
                ncol(bulk_counts), nrow(bulk_counts))
    log_message("样本分组分布:")
    print(table(sample_info$group))
    log_message("批次分布:")
    print(table(sample_info$batch))
    min_size <- min(sample_info$library_size)
    max_size <- max(sample_info$library_size)
    log_message("文库大小范围: %d - %d", min_size, max_size)
  }
  
  # 准备基因元数据
  gene_metadata <- data.frame(
    gene_id = rownames(bulk_counts)
  )
  
  # 如果单细胞数据中有行数据，则合并
  if (!is.null(rowData(sce)) && ncol(rowData(sce)) > 0) {
    gene_info <- as.data.frame(rowData(sce))
    gene_metadata <- cbind(gene_metadata, gene_info[, -1, drop = FALSE])
  }
  
  # 返回结果
  result <- list(
    bulk_counts = bulk_counts,
    bulk_metadata = sample_info,
    gene_metadata = gene_metadata,
    cell_proportions = cell_proportions,
    simulation_params = list(
      confounding_level = confounding_level,
      batch_effect_strength = batch_effect_strength,
      n_batches = n_batches
    )
  )
  
  return(result)
}

# 模块5: 主控制器函数 ------------------------------------------------------
#' @description 运行整个模拟管道，生成多个模拟数据集
#' @param param_grid 包含所有参数组合的数据框
#' @param output_dir 输出目录
#' @param n_cores 用于并行计算的核心数
#' @param common_params 所有模拟共享的参数
#' @param save_intermediate 是否保存中间单细胞数据
#' @param verbose 是否打印详细信息
#' @return 生成的所有模拟数据集的列表
run_simulation_pipeline <- function(
    param_grid = NULL,
    output_dir = "simulation_results",
    n_cores = 1,
    common_params = list(
      n_cells = 3000,
      n_genes = 3000,
      n_groups = 6,
      n_bulk_samples = 30
    ),
    save_intermediate = FALSE,
    verbose = TRUE
) {
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 如果没有提供参数网格，则创建默认的参数网格
  if (is.null(param_grid)) {
    param_grid <- expand.grid(
      intra_group_distance = c(0.1, 0.3, 0.5),       # 近、中、远
      inter_group_distance = c(0.5, 1.5, 3.0),       # 近、中、远
      overall_cell_spread = c(0.1, 0.3, 0.5),        # 近、中、远
      group_separation_strength = c(0.8, 1.2, 1.6),  # 弱、中、强
      confounding_level = c("low", "medium", "high"),
      batch_effect_strength = c(0.3, 0.6, 0.9),      # 弱、中、强
      seed = 42:(42 + n_cores - 1),                  # 确保不同核心使用不同种子
      stringsAsFactors = FALSE
    )
  }
  
  if (verbose) {
    log_message("将生成%d个不同参数组合的模拟数据集", nrow(param_grid))
    log_message("使用%d个核心进行并行计算", n_cores)
  }
  
  # 设置并行后端
  if (n_cores > 1) {
    bp_param <- MulticoreParam(workers = n_cores)
  } else {
    bp_param <- SerialParam()
  }
  
  # 定义单次模拟的函数
  run_single_simulation <- function(params_row, row_idx) {
    # 提取这次迭代的参数
    iter_params <- as.list(params_row)
    iter_params <- c(iter_params, common_params)
    
    # 创建唯一的文件名前缀
    prefix <- sprintf(
      "sim_intra%.1f_inter%.1f_spread%.1f_sep%.1f_conf_%s_batch%.1f_seed%d",
      iter_params$intra_group_distance,
      iter_params$inter_group_distance,
      iter_params$overall_cell_spread,
      iter_params$group_separation_strength,
      iter_params$confounding_level,
      iter_params$batch_effect_strength,
      iter_params$seed
    )
    
    # 替换文件名中的点
    prefix <- gsub("\\.", "p", prefix)
    
    if (verbose) {
      log_message("开始第%d/%d次模拟: %s", 
                  row_idx, nrow(param_grid), prefix)
    }
    
    tryCatch({
      # 步骤1: 配置参数
      sim_params <- configure_simulation_params(
        intra_group_distance = iter_params$intra_group_distance,
        inter_group_distance = iter_params$inter_group_distance,
        overall_cell_spread = iter_params$overall_cell_spread,
        group_separation_strength = iter_params$group_separation_strength,
        n_cells = iter_params$n_cells,
        n_genes = iter_params$n_genes,
        n_groups = iter_params$n_groups,
        seed = iter_params$seed,
        verbose = verbose
      )
      
      # 步骤1.5: 先生成临时数据以获取基因名称格式
      temp_sc_data <- splatSimulate(
        sim_params$splatter_params,
        method = "groups",
        verbose = FALSE
      )
      actual_gene_names <- rownames(temp_sc_data)
      
      # 步骤2: 生成单细胞数据
      # 创建几个共表达模块作为ground truth
      n_modules <- 3
      n_genes_per_module <- 50
      
      # 使用实际基因名创建共表达模块
      co_expression_modules <- list()
      for (i in 1:n_modules) {
        selected_genes <- sample(actual_gene_names, n_genes_per_module)
        co_expression_modules[[paste0("module", i)]] <- selected_genes
      }
      
      # 生成单细胞数据
      sc_data <- generate_single_cell_data(
        params = sim_params,
        co_expression_modules = co_expression_modules,
        seed = iter_params$seed,
        verbose = verbose
      )
      
      # 步骤3: 增强共表达网络
      # 计算共表达强度，基于组间距离和组内距离
      co_expr_strength <- 0.5 + 
        (iter_params$intra_group_distance * 0.5) - 
        (iter_params$inter_group_distance * 0.1)
      co_expr_strength <- max(0.1, min(0.9, co_expr_strength))
      
      enhanced_sc_data <- enhance_co_expression(
        sce = sc_data,
        celltype_col = "CellType",
        assay_use = "counts",
        co_expression_strength = co_expr_strength,
        verbose = verbose
      )
      
      # 步骤4: 聚合为bulk样本
      bulk_data <- aggregate_to_bulk(
        sce = enhanced_sc_data,
        n_bulk_samples = iter_params$n_bulk_samples,
        confounding_level = iter_params$confounding_level,
        batch_effect_strength = iter_params$batch_effect_strength,
        seed = iter_params$seed,
        verbose = verbose
      )
      
      # 保存结果
      output_file <- file.path(output_dir, paste0(prefix, ".rds"))
      saveRDS(bulk_data, output_file)
      
      # 可选：保存中间的单细胞数据
      if (save_intermediate) {
        sc_output_file <- file.path(output_dir, paste0(prefix, "_sc_data.rds"))
        saveRDS(enhanced_sc_data, sc_output_file)
      }
      
      if (verbose) {
        log_message("第%d/%d次模拟完成: %s", 
                    row_idx, nrow(param_grid), prefix)
      }
      
      # 返回成功状态
      return(list(
        status = "success",
        file = output_file,
        params = params_row
      ))
      
    }, error = function(e) {
      log_message("第%d/%d次模拟失败: %s", 
                  row_idx, nrow(param_grid), as.character(e),
                  level = "ERROR")
      
      # 返回错误状态
      return(list(
        status = "error",
        error_message = as.character(e),
        params = params_row
      ))
    })
  }
  
  # 并行运行所有模拟
  start_time <- Sys.time()
  
  results <- bplapply(
    X = seq_len(nrow(param_grid)),
    FUN = function(i) {
      run_single_simulation(param_grid[i, ], i)
    },
    BPPARAM = bp_param
  )
  
  end_time <- Sys.time()
  
  # 汇总结果
  success_count <- sum(sapply(results, function(x) x$status == "success"))
  error_count <- nrow(param_grid) - success_count
  
  if (verbose) {
    log_message("模拟完成. 成功: %d, 失败: %d", success_count, error_count)
    log_message("总运行时间: %s", format(end_time - start_time))
  }
  
  # 保存运行概要
  summary_file <- file.path(output_dir, "simulation_summary.rds")
  summary_data <- list(
    param_grid = param_grid,
    results = results,
    runtime = end_time - start_time,
    success_count = success_count,
    error_count = error_count,
    timestamp = Sys.time()
  )
  saveRDS(summary_data, summary_file)
  
  return(invisible(results))
}

# 模块6: 辅助功能和可视化函数 ---------------------------------------------
#' @description 可视化bulk RNA-seq数据的批次效应和组效应 (修复版)
#' @param bulk_data 由aggregate_to_bulk生成的bulk数据
#' @param plot_file 输出文件名，如果为NULL则不保存
#' @param filter_zero_var 是否过滤零方差基因
#' @return 包含多个图的列表
visualize_bulk_data <- function(bulk_data, plot_file = NULL, filter_zero_var = TRUE) {
  # 提取数据
  counts <- bulk_data$bulk_counts
  metadata <- bulk_data$bulk_metadata
  
  # 确保行名和列名正确设置
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("Gene", 1:nrow(counts))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("Sample", 1:ncol(counts))
  }
  
  # 执行VST转换
  vsd <- NULL
  tryCatch({
    message("converting counts to integer mode")
    vsd <- vst(round(counts))
  }, error = function(e) {
    # 如果VST失败，尝试使用简单的log转换
    warning("VST转换失败，使用简单的log转换代替: ", e$message)
    pseudo_count <- 1
    vsd <- log2(counts + pseudo_count)
  })
  
  # 过滤零方差基因以避免PCA错误
  if (filter_zero_var) {
    gene_vars <- apply(vsd, 1, var)
    zero_var_genes <- which(gene_vars < 1e-10)  # 接近零的方差
    
    if (length(zero_var_genes) > 0) {
      message(sprintf("过滤了%d个零/低方差基因以避免PCA错误", length(zero_var_genes)))
      if (length(zero_var_genes) < nrow(vsd)) {
        vsd <- vsd[-zero_var_genes, ]
      } else {
        # 如果所有基因都有零方差，添加微小的噪声
        warning("所有基因都有近零方差，添加微小噪声以执行PCA")
        vsd <- vsd + matrix(rnorm(length(vsd), 0, 0.01), nrow = nrow(vsd))
      }
    }
  }
  
  # 执行PCA
  pca_data <- prcomp(t(vsd), scale. = TRUE)
  
  # 创建可视化数据框
  pca_df <- data.frame(
    PC1 = pca_data$x[, 1],
    PC2 = pca_data$x[, 2],
    stringsAsFactors = FALSE
  )
  
  # 检查metadata中的列名 - 解决大小写问题
  group_col <- NULL
  batch_col <- NULL
  
  # 检测组和批次列的名称（考虑大小写）
  if ("group" %in% colnames(metadata)) {
    group_col <- "group"
  } else if ("Group" %in% colnames(metadata)) {
    group_col <- "Group"
  } else {
    group_col <- "group"
    metadata$group <- rep("Unknown", nrow(metadata))
    warning("未找到group或Group列，添加默认值")
  }
  
  if ("batch" %in% colnames(metadata)) {
    batch_col <- "batch"
  } else if ("Batch" %in% colnames(metadata)) {
    batch_col <- "Batch"
  } else {
    batch_col <- "batch"
    metadata$batch <- rep("Unknown", nrow(metadata))
    warning("未找到batch或Batch列，添加默认值")
  }
  
  # 添加组和批次到PCA数据框
  pca_df$Group <- metadata[[group_col]]
  pca_df$Batch <- metadata[[batch_col]]
  
  # 创建PCA图 - 按组着色
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(color = Group), level = 0.95) +
    theme_bw() +
    labs(title = "PCA: 按组着色")
  
  # 创建PCA图 - 按批次着色
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(color = Batch), level = 0.95) +
    theme_bw() +
    labs(title = "PCA: 按批次着色")
  
  # 创建文库大小分布图 - 使用检测到的列名
  p3 <- ggplot(metadata, aes_string(x = paste0("factor(", batch_col, ")"), 
                                    y = "library_size", 
                                    fill = paste0("factor(", group_col, ")"))) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = "文库大小分布", 
         y = "文库大小",
         x = "批次",
         fill = "组别") 
  
  # 创建细胞类型比例图
  cell_types <- grep("^prop_", colnames(metadata), value = TRUE)
  
  # 确保有细胞类型比例信息
  if (length(cell_types) > 0) {
    prop_data <- reshape2::melt(metadata, 
                                id.vars = c("sample_id", group_col, batch_col),
                                measure.vars = cell_types,
                                variable.name = "cell_type",
                                value.name = "proportion")
    prop_data$cell_type <- gsub("prop_", "", prop_data$cell_type)
    
    # 使用检测到的组列名
    p4 <- ggplot(prop_data, aes_string(x = group_col, y = "proportion", fill = "cell_type")) +
      geom_boxplot(position = "dodge") +
      theme_bw() +
      labs(title = "各组中的细胞类型比例", 
           y = "比例", 
           x = "组别",
           fill = "细胞类型")
  } else {
    # 如果没有细胞类型数据，创建一个空白图
    p4 <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, 
               label = "没有找到细胞类型比例数据\n(prop_* 列)", size = 5)
  }
  
  # 组合所有图
  plots <- list(
    pca_by_group = p1,
    pca_by_batch = p2,
    library_size = p3,
    cell_proportions = p4
  )
  
  # 如果提供了文件名，则保存图
  if (!is.null(plot_file)) {
    pdf(plot_file, width = 12, height = 10)
    gridExtra::grid.arrange(
      p1, p2, p3, p4, 
      ncol = 2,
      top = "Bulk RNA-seq 数据质量控制"
    )
    dev.off()
    message(sprintf("图表已保存到: %s", plot_file))
  }
  
  return(plots)
}

#' @description 使用常见的批次校正方法校正数据并评估校正效果 (简化版)
#' @param bulk_data 由aggregate_to_bulk生成的bulk数据
#' @param methods 要使用的方法向量 ("ComBat" 和/或 "limma")
#' @param plot_file 输出文件名，如果为NULL则不保存
#' @param filter_zero_var 是否过滤零方差基因
#' @return 校正后的数据列表
batch_correction_benchmark <- function(
    bulk_data,
    methods = c("ComBat", "limma"),  # 移除了"RUV"选项
    plot_file = NULL,
    filter_zero_var = TRUE
) {
  # 加载必要的包
  suppressPackageStartupMessages({
    require(sva)       # for ComBat
    require(limma)     # for removeBatchEffect
  })
  
  # 提取数据
  counts <- bulk_data$bulk_counts
  metadata <- bulk_data$bulk_metadata
  
  # 确保行名和列名正确设置
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("Gene", 1:nrow(counts))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("Sample", 1:ncol(counts))
  }
  
  # 检查metadata中的列名 - 解决大小写问题
  group_col <- NULL
  batch_col <- NULL
  
  # 检测组和批次列的名称（考虑大小写）
  if ("group" %in% colnames(metadata)) {
    group_col <- "group"
  } else if ("Group" %in% colnames(metadata)) {
    group_col <- "Group"
  } else {
    group_col <- "group"
    metadata$group <- rep("Unknown", nrow(metadata))
    warning("未找到group或Group列，添加默认值")
  }
  
  if ("batch" %in% colnames(metadata)) {
    batch_col <- "batch"
  } else if ("Batch" %in% colnames(metadata)) {
    batch_col <- "Batch"
  } else {
    batch_col <- "batch"
    metadata$batch <- rep("Unknown", nrow(metadata))
    warning("未找到batch或Batch列，添加默认值")
  }
  
  # 执行VST转换
  vsd <- NULL
  tryCatch({
    message("converting counts to integer mode")
    vsd <- vst(round(counts))
  }, error = function(e) {
    # 如果VST失败，尝试使用简单的log转换
    warning("VST转换失败，使用简单的log转换代替: ", e$message)
    pseudo_count <- 1
    vsd <- log2(counts + pseudo_count)
  })
  
  # 过滤零方差基因以避免PCA错误
  if (filter_zero_var) {
    gene_vars <- apply(vsd, 1, var)
    zero_var_genes <- which(gene_vars < 1e-10)  # 接近零的方差
    
    if (length(zero_var_genes) > 0) {
      message(sprintf("过滤了%d个零/低方差基因以避免分析错误", length(zero_var_genes)))
      if (length(zero_var_genes) < nrow(vsd)) {
        vsd <- vsd[-zero_var_genes, ]
      } else {
        # 如果所有基因都有零方差，添加微小的噪声
        warning("所有基因都有近零方差，添加微小噪声以执行分析")
        vsd <- vsd + matrix(rnorm(length(vsd), 0, 0.01), nrow = nrow(vsd))
      }
    }
  }
  
  # 创建存储校正后数据的列表
  corrected_data <- list(
    original = vsd
  )
  
  # 应用批次校正方法
  if ("ComBat" %in% methods) {
    # ComBat校正
    tryCatch({
      combat_data <- ComBat(
        dat = as.matrix(vsd),
        batch = metadata[[batch_col]],
        mod = model.matrix(as.formula(paste0("~", group_col)), data = metadata)
      )
      corrected_data$ComBat <- combat_data
    }, error = function(e) {
      warning("ComBat校正失败: ", e$message)
    })
  }
  
  if ("limma" %in% methods) {
    # limma校正
    tryCatch({
      limma_data <- removeBatchEffect(
        x = as.matrix(vsd),
        batch = metadata[[batch_col]],
        design = model.matrix(as.formula(paste0("~", group_col)), data = metadata)
      )
      corrected_data$limma <- limma_data
    }, error = function(e) {
      warning("limma校正失败: ", e$message)
    })
  }
  
  # 添加简单的均值中心化批次校正
  tryCatch({
    message("使用简单的批次均值中心化方法...")
    
    # 获取批次
    batches <- unique(metadata[[batch_col]])
    
    # 对每个批次进行均值中心化
    simple_corrected <- vsd
    for (b in batches) {
      batch_samples <- which(metadata[[batch_col]] == b)
      batch_means <- rowMeans(vsd[, batch_samples, drop=FALSE])
      # 减去批次均值，加上总体均值
      simple_corrected[, batch_samples] <- vsd[, batch_samples, drop=FALSE] - batch_means + rowMeans(vsd)
    }
    
    corrected_data$SimpleBatchCorrection <- simple_corrected
  }, error = function(e) {
    warning("简单批次校正失败: ", e$message)
  })
  
  # 对每种校正方法进行PCA分析并可视化
  if (!is.null(plot_file)) {
    # 检查有多少成功的方法
    success_methods <- names(corrected_data)
    n_methods <- length(success_methods)
    
    if (n_methods > 0) {
      pdf(plot_file, width = 15, height = 10)
      par(mfrow = c(2, ceiling(n_methods / 2)))
      
      for (method_name in success_methods) {
        # 执行PCA
        method_data <- corrected_data[[method_name]]
        
        # 过滤零方差列
        if (filter_zero_var) {
          col_vars <- apply(t(method_data), 2, var)
          zero_var_cols <- col_vars < 1e-10
          if (any(zero_var_cols)) {
            if (!all(zero_var_cols)) {
              method_data <- method_data[!zero_var_cols, ]
            } else {
              # 如果所有列都有零方差，添加微小的噪声
              method_data <- method_data + matrix(rnorm(length(method_data), 0, 0.01), 
                                                  nrow = nrow(method_data))
            }
          }
        }
        
        # 尝试执行PCA
        tryCatch({
          pca_data <- prcomp(t(method_data), scale. = TRUE)
          
          # 创建可视化数据框
          pca_df <- data.frame(
            PC1 = pca_data$x[, 1],
            PC2 = pca_data$x[, 2],
            Group = metadata[[group_col]],
            Batch = metadata[[batch_col]]
          )
          
          # 绘制PCA图
          plot(
            pca_df$PC1, pca_df$PC2,
            col = as.numeric(factor(pca_df$Batch)),
            pch = as.numeric(factor(pca_df$Group)),
            main = paste("PCA -", method_name),
            xlab = "PC1", ylab = "PC2"
          )
          
          # 添加图例
          legend("topright", 
                 legend = c(
                   levels(factor(pca_df$Group)), 
                   levels(factor(pca_df$Batch))
                 ),
                 pch = c(rep(1:length(unique(pca_df$Group)), 1),
                         rep(1, length(unique(pca_df$Batch)))),
                 col = c(rep(1, length(unique(pca_df$Group))), 
                         1:length(unique(pca_df$Batch))),
                 title = "Group & Batch")
        }, error = function(e) {
          # 如果PCA失败，创建一个空白图
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, paste("PCA failed for", method_name, ":", e$message), cex = 0.8)
        })
      }
      
      dev.off()
      message(sprintf("批次校正比较图已保存到: %s", plot_file))
    } else {
      warning("没有成功的校正方法可视化")
    }
  }
  
  return(corrected_data)
}
# 提供给用户的主要函数 ----------------------------------------------------
#' @description 运行RNA-seq模拟主程序 (修复scan模式错误)
#' @param mode 运行模式: "simple", "scan", "custom"
#' @param params 自定义参数列表
#' @param output_dir 输出目录
#' @param n_cores 用于并行计算的核心数
#' @param verbose 是否打印详细信息
#' @return 生成的模拟数据或运行结果
run_simulation <- function(
    mode = "simple",
    params = NULL,
    output_dir = "simulation_results",
    n_cores = 1,
    verbose = TRUE
) {
  if (mode == "simple") {
    # 简单模式: 运行单个示例模拟
    if (verbose) log_message("运行简单示例模拟...")
    
    # 配置参数
    sim_params <- configure_simulation_params(
      intra_group_distance = 0.2,
      inter_group_distance = 2.0,
      overall_cell_spread = 0.3,
      group_separation_strength = 1.2,
      verbose = verbose
    )
    
    # 重要修复：先生成单细胞数据以获取正确的基因名格式
    temp_sc_data <- splatSimulate(
      sim_params$splatter_params,
      method = "groups",
      verbose = FALSE
    )
    
    # 获取实际基因名
    actual_gene_names <- rownames(temp_sc_data)
    if (verbose) {
      log_message("检测到实际基因名称格式: %s, %s, ...", 
                  actual_gene_names[1], actual_gene_names[2])
    }
    
    # 创建共表达模块（使用实际基因名）
    set.seed(42)
    n_modules <- 3
    n_genes_per_module <- 50
    
    co_expression_modules <- list()
    for (i in 1:n_modules) {
      # 直接从实际基因名中随机选择
      selected_genes <- sample(actual_gene_names, n_genes_per_module)
      co_expression_modules[[paste0("module", i)]] <- selected_genes
      
      if (verbose) {
        log_message("为模块%d选择了%d个基因: %s, %s, ...", 
                    i, length(selected_genes), selected_genes[1], selected_genes[2])
      }
    }
    
    # 生成单细胞数据（使用正确格式的基因模块）
    sc_data <- generate_single_cell_data(
      params = sim_params,
      co_expression_modules = co_expression_modules,
      seed = 42,
      verbose = verbose
    )
    
    # 检查模块分配是否成功
    module_cols <- grep("^in_module_", colnames(rowData(sc_data)), value = TRUE)
    if (length(module_cols) > 0) {
      for (col in module_cols) {
        n_genes_in_module <- sum(rowData(sc_data)[[col]])
        if (verbose) {
          log_message("模块 %s 包含 %d 个基因", 
                      gsub("^in_module_", "", col), n_genes_in_module)
        }
      }
    }
    
    # 增强共表达网络
    enhanced_sc_data <- enhance_co_expression(
      sce = sc_data,
      celltype_col = "CellType",
      assay_use = "counts",
      co_expression_strength = 0.7,
      verbose = verbose
    )
    
    # 聚合为bulk样本
    bulk_data <- aggregate_to_bulk(
      sce = enhanced_sc_data,
      n_bulk_samples = 30,
      confounding_level = "medium",
      batch_effect_strength = 0.6,
      seed = 42,
      verbose = verbose
    )
    
    # 最终检查生成的数据中的模块
    module_cols <- grep("^in_module_", colnames(bulk_data$gene_metadata), value = TRUE)
    if (length(module_cols) > 0) {
      for (col in module_cols) {
        n_genes_in_module <- sum(bulk_data$gene_metadata[[col]])
        if (verbose) {
          log_message("最终结果: 模块 %s 包含 %d 个基因", 
                      gsub("^in_module_", "", col), n_genes_in_module)
        }
      }
    }
    
    # 保存结果
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, "simple_example.rds")
    saveRDS(bulk_data, output_file)
    
    result <- list(
      bulk_data = bulk_data,
      output_file = output_file
    )
    
    return(result)
  } else if (mode == "scan") {
    # 扫描模式: 运行参数扫描
    if (verbose) log_message("运行参数扫描模拟...")
    
    # 创建小型参数网格进行扫描
    mini_grid <- expand.grid(
      intra_group_distance = c(0.1, 0.3),
      inter_group_distance = c(1.0, 3.0),
      overall_cell_spread = c(0.2),
      group_separation_strength = c(1.0),
      confounding_level = c("low", "high"),
      batch_effect_strength = c(0.5),
      seed = 42:(42 + n_cores - 1),
      stringsAsFactors = FALSE
    )
    
    # 运行管道
    scan_result <- run_simulation_pipeline(  # 修复：存储结果到 scan_result
      param_grid = mini_grid,
      output_dir = output_dir,
      n_cores = n_cores,
      common_params = list(
        n_cells = 1000,  # 减少细胞数以加快测试
        n_genes = 1000,  # 减少基因数以加快测试
        n_groups = 4,    # 减少组数以简化测试
        n_bulk_samples = 20
      ),
      verbose = verbose
    )
    
    return(scan_result)  # 修复：返回 scan_result
  } else if (mode == "custom") {
    # 自定义模式: 使用用户提供的参数
    if (is.null(params)) {
      stop("自定义模式需要提供params参数")
    }
    
    if (verbose) log_message("运行自定义参数模拟...")
    
    # 检查参数结构并运行相应的函数
    if ("param_grid" %in% names(params)) {
      # 如果提供了参数网格，则运行参数扫描
      result <- run_simulation_pipeline(
        param_grid = params$param_grid,
        output_dir = output_dir,
        n_cores = n_cores,
        common_params = params$common_params,
        verbose = verbose
      )
    } else {
      # 否则运行单个模拟
      result <- list()
      
      # 修复: 分离参数
      # 提取configure_simulation_params需要的参数
      config_params_names <- names(formals(configure_simulation_params))
      config_params <- params[names(params) %in% config_params_names]
      
      # 保存不属于configure_simulation_params的参数
      other_params <- params[!names(params) %in% config_params_names]
      
      # 配置参数
      sim_params <- do.call(configure_simulation_params, config_params)
      
      # 重要：先生成临时数据以获取正确的基因名格式
      temp_sc_data <- splatSimulate(
        sim_params$splatter_params,
        method = "groups",
        verbose = FALSE
      )
      actual_gene_names <- rownames(temp_sc_data)
      
      # 创建共表达模块（使用实际基因名）
      set.seed(params$seed)
      n_modules <- 3
      n_genes_per_module <- 50
      
      co_expression_modules <- list()
      for (i in 1:n_modules) {
        selected_genes <- sample(actual_gene_names, n_genes_per_module)
        co_expression_modules[[paste0("module", i)]] <- selected_genes
      }
      
      # 如果用户提供了自定义模块，则使用用户提供的模块
      if (!is.null(params$co_expression_modules)) {
        co_expression_modules <- params$co_expression_modules
      }
      
      # 生成单细胞数据
      sc_data <- generate_single_cell_data(
        params = sim_params,
        co_expression_modules = co_expression_modules,
        seed = params$seed,
        verbose = verbose
      )
      
      # 增强共表达网络
      co_expression_strength <- 0.7  # 默认值
      if (!is.null(other_params$co_expression_strength)) {
        co_expression_strength <- other_params$co_expression_strength
      }
      
      enhanced_data <- enhance_co_expression(
        sce = sc_data,
        co_expression_strength = co_expression_strength,
        verbose = verbose
      )
      
      # 聚合为bulk样本
      bulk_params <- list(
        sce = enhanced_data,
        seed = params$seed,
        verbose = verbose
      )
      
      # 添加其他适用于aggregate_to_bulk的参数
      if (!is.null(other_params$n_bulk_samples)) 
        bulk_params$n_bulk_samples <- other_params$n_bulk_samples
      
      if (!is.null(other_params$confounding_level)) 
        bulk_params$confounding_level <- other_params$confounding_level
      
      if (!is.null(other_params$batch_effect_strength)) 
        bulk_params$batch_effect_strength <- other_params$batch_effect_strength
      
      if (!is.null(other_params$n_batches)) 
        bulk_params$n_batches <- other_params$n_batches
      
      # 调用函数时使用do.call以正确处理参数
      bulk_data <- do.call(aggregate_to_bulk, bulk_params)
      
      result$bulk_data <- bulk_data
      
      # 保存结果
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- file.path(output_dir, "custom_simulation.rds")
        saveRDS(bulk_data, output_file)
        result$output_file <- output_file
      }
    }
    
    return(result)
  } else {
    stop("无效的模式。支持的模式有: 'simple', 'scan', 'custom'")
  }
}
