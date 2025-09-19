#' @title 高质量Pseudo-Bulk RNA-seq数据模拟管道 (黄金标准版)
#' @description 基于黄金标准原则生成具有生物学和统计学真实性的pseudo-bulk数据
#' @author 修订版 - 专攻转录组数据分析
#' @date 2025-01-19
#' 
#' 本脚本实现了以下黄金标准特征:
#' 1. 生物学真实性: 细胞类型特异性共表达模块
#' 2. 统计学真实性: 多样化抽样策略和可变测序深度
#' 3. 参数化设计: 高度可配置的模拟参数
#' 4. 模块化架构: 便于扩展和维护

# 必要库的加载 -------------------------------------------------------------
suppressPackageStartupMessages({
  library(BiocParallel)       # 并行计算
  library(splatter)           # 单细胞RNA-seq模拟
  library(scDesign3)          # 共表达网络增强
  library(SingleCellExperiment) # 单细胞数据处理
  library(scater)             # 单细胞数据分析
  library(Matrix)             # 稀疏矩阵处理
  library(ggplot2)            # 可视化
  library(dplyr)              # 数据操作
  library(WGCNA)              # 基因共表达网络分析
  library(igraph)             # 网络分析
})

# 辅助函数: 增强的日志系统 -------------------------------------------------
log_message <- function(msg, level = "INFO", ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (...length() > 0) {
    msg <- sprintf(msg, ...)
  }
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

# 模块1: 细胞类型特异性共表达模块设计器 ------------------------------------
#' @description 设计细胞类型特异性的共表达模块
#' @param n_cell_types 细胞类型数量
#' @param n_genes_total 总基因数量
#' @param modules_per_celltype 每个细胞类型的模块数量
#' @param genes_per_module 每个模块的基因数量范围
#' @param shared_modules_fraction 跨细胞类型共享模块的比例
#' @param correlation_strength 模块内基因相关性强度
#' @param seed 随机种子
#' @return 包含细胞类型特异性模块信息的列表
design_celltype_specific_modules <- function(
    n_cell_types = 6,
    n_genes_total = 3000,
    modules_per_celltype = 2,
    genes_per_module = c(30, 80),  # 最小和最大基因数
    shared_modules_fraction = 0.2,  # 20%的模块是跨细胞类型共享的
    correlation_strength = c(0.6, 0.9),  # 相关性强度范围
    seed = 42
) {
  set.seed(seed)
  
  log_message("设计细胞类型特异性共表达模块...")
  
  # 生成基因名称
  gene_names <- paste0("Gene", sprintf("%04d", 1:n_genes_total))
  
  # 初始化模块结构
  celltype_modules <- list()
  shared_modules <- list()
  gene_assignments <- data.frame(
    gene_id = gene_names,
    stringsAsFactors = FALSE
  )
  
  # 为每个细胞类型添加列
  for (i in 1:n_cell_types) {
    celltype_name <- paste0("CellType", i)
    gene_assignments[[paste0("in_", celltype_name)]] <- FALSE
    celltype_modules[[celltype_name]] <- list()
  }
  
  # 计算共享模块数量
  total_modules <- n_cell_types * modules_per_celltype
  n_shared_modules <- round(total_modules * shared_modules_fraction)
  n_specific_modules <- total_modules - n_shared_modules
  
  used_genes <- c()  # 跟踪已使用的基因
  
  # 第一步: 创建共享模块（跨多个细胞类型）
  if (n_shared_modules > 0) {
    log_message("创建%d个跨细胞类型共享模块...", n_shared_modules)
    
    for (i in 1:n_shared_modules) {
      # 随机选择要共享的细胞类型（2-4个）
      n_sharing_types <- sample(2:min(4, n_cell_types), 1)
      sharing_types <- sample(1:n_cell_types, n_sharing_types)
      
      # 确定模块大小
      module_size <- sample(genes_per_module[1]:genes_per_module[2], 1)
      
      # 选择基因（避免重复使用）
      available_genes <- setdiff(gene_names, used_genes)
      if (length(available_genes) < module_size) {
        warning("可用基因不足，减少模块大小")
        module_size <- length(available_genes)
      }
      
      module_genes <- sample(available_genes, module_size)
      used_genes <- c(used_genes, module_genes)
      
      # 为模块分配相关性强度
      corr_strength <- runif(1, correlation_strength[1], correlation_strength[2])
      
      # 存储共享模块信息
      module_name <- paste0("SharedModule", i)
      shared_modules[[module_name]] <- list(
        genes = module_genes,
        cell_types = paste0("CellType", sharing_types),
        correlation_strength = corr_strength,
        module_type = "shared"
      )
      
      # 更新基因分配表
      for (ct_idx in sharing_types) {
        celltype_name <- paste0("CellType", ct_idx)
        gene_assignments[gene_assignments$gene_id %in% module_genes, 
                         paste0("in_", celltype_name)] <- TRUE
        
        # 添加到细胞类型特异模块列表
        celltype_modules[[celltype_name]][[module_name]] <- list(
          genes = module_genes,
          correlation_strength = corr_strength,
          module_type = "shared"
        )
      }
    }
  }
  
  # 第二步: 创建细胞类型特异性模块
  log_message("创建细胞类型特异性模块...")
  
  specific_modules_per_type <- n_specific_modules %/% n_cell_types
  remaining_modules <- n_specific_modules %% n_cell_types
  
  for (i in 1:n_cell_types) {
    celltype_name <- paste0("CellType", i)
    
    # 确定此细胞类型的特异性模块数量
    n_modules_this_type <- specific_modules_per_type
    if (i <= remaining_modules) {
      n_modules_this_type <- n_modules_this_type + 1
    }
    
    for (j in 1:n_modules_this_type) {
      # 确定模块大小
      module_size <- sample(genes_per_module[1]:genes_per_module[2], 1)
      
      # 选择基因
      available_genes <- setdiff(gene_names, used_genes)
      if (length(available_genes) < module_size) {
        if (length(available_genes) == 0) break
        module_size <- length(available_genes)
      }
      
      module_genes <- sample(available_genes, module_size)
      used_genes <- c(used_genes, module_genes)
      
      # 分配相关性强度（特异性模块通常更强）
      corr_strength <- runif(1, 
                             correlation_strength[1] + 0.1, 
                             correlation_strength[2])
      
      # 存储模块信息
      module_name <- paste0(celltype_name, "_SpecificModule", j)
      celltype_modules[[celltype_name]][[module_name]] <- list(
        genes = module_genes,
        correlation_strength = corr_strength,
        module_type = "specific"
      )
      
      # 更新基因分配表
      gene_assignments[gene_assignments$gene_id %in% module_genes, 
                       paste0("in_", celltype_name)] <- TRUE
    }
  }
  
  # 统计信息
  total_assigned_genes <- length(used_genes)
  unassigned_genes <- n_genes_total - total_assigned_genes
  
  log_message("模块设计完成:")
  log_message("  - 总模块数: %d (共享: %d, 特异: %d)", 
              length(shared_modules) + sum(sapply(celltype_modules, length)) - length(shared_modules),
              length(shared_modules), 
              sum(sapply(celltype_modules, length)) - length(shared_modules))
  log_message("  - 已分配基因: %d/%d (%.1f%%)", 
              total_assigned_genes, n_genes_total, 
              100 * total_assigned_genes / n_genes_total)
  
  # 返回结果
  return(list(
    celltype_modules = celltype_modules,
    shared_modules = shared_modules,
    gene_assignments = gene_assignments,
    design_params = list(
      n_cell_types = n_cell_types,
      n_genes_total = n_genes_total,
      modules_per_celltype = modules_per_celltype,
      genes_per_module = genes_per_module,
      shared_modules_fraction = shared_modules_fraction,
      correlation_strength = correlation_strength
    )
  ))
}

# 模块2: 增强的单细胞数据生成器 ----------------------------------------------
#' @description 生成具有预定义共表达结构的单细胞数据
#' @param module_design 由design_celltype_specific_modules返回的模块设计
#' @param n_cells 细胞总数
#' @param group_probabilities 各细胞类型比例
#' @param de_gene_proportion 差异表达基因比例
#' @param seed 随机种子
#' @param verbose 是否打印详细信息
#' @return 增强的SingleCellExperiment对象
generate_enhanced_single_cell_data <- function(
    module_design,
    n_cells = 3000,
    group_probabilities = NULL,
    de_gene_proportion = 0.2,
    seed = 42,
    verbose = TRUE
) {
  set.seed(seed)
  
  if (verbose) log_message("生成具有预定义共表达结构的单细胞数据...")
  
  n_cell_types <- module_design$design_params$n_cell_types
  n_genes <- module_design$design_params$n_genes_total
  
  # 设置细胞类型比例
  if (is.null(group_probabilities)) {
    group_probabilities <- rep(1/n_cell_types, n_cell_types)
  }
  
  # 配置Splatter参数
  params <- newSplatParams(
    batchCells = n_cells,
    nGenes = n_genes,
    group.prob = group_probabilities,
    de.prob = rep(de_gene_proportion, n_cell_types),
    de.facLoc = 1.0,
    seed = seed
  )
  
  # 生成基础单细胞数据
  if (verbose) log_message("使用Splatter生成基础单细胞表达谱...")
  sim <- splatSimulate(params, method = "groups", verbose = FALSE)
  
  # 设置细胞类型标签
  cell_type_levels <- paste0("CellType", 1:n_cell_types)
  if (all(grepl("^Group", sim$Group))) {
    colData(sim)$CellType <- factor(
      gsub("Group", "CellType", sim$Group), 
      levels = cell_type_levels
    )
  } else {
    colData(sim)$CellType <- factor(
      sim$Group, 
      levels = 1:n_cell_types, 
      labels = cell_type_levels
    )
  }
  
  # 应用细胞类型特异性共表达结构
  if (verbose) log_message("应用细胞类型特异性共表达结构...")
  
  enhanced_counts <- counts(sim)
  
  # 对每个细胞类型分别处理
  for (celltype_name in names(module_design$celltype_modules)) {
    if (verbose) log_message("处理%s的共表达模块...", celltype_name)
    
    # 获取该细胞类型的细胞
    celltype_cells <- which(sim$CellType == celltype_name)
    if (length(celltype_cells) == 0) next
    
    # 获取该细胞类型的模块
    celltype_modules <- module_design$celltype_modules[[celltype_name]]
    
    # 对每个模块应用共表达结构
    for (module_name in names(celltype_modules)) {
      module_info <- celltype_modules[[module_name]]
      module_genes <- module_info$genes
      corr_strength <- module_info$correlation_strength
      
      # 检查基因是否存在于数据中
      valid_genes <- intersect(module_genes, rownames(enhanced_counts))
      if (length(valid_genes) < 2) next
      
      # 获取当前表达水平
      current_expr <- enhanced_counts[valid_genes, celltype_cells, drop = FALSE]
      
      # 应用共表达结构使用改进的方法
      enhanced_expr <- apply_coexpression_structure(
        expr_matrix = current_expr,
        correlation_strength = corr_strength,
        method = "factor_model"
      )
      
      # 更新表达矩阵
      enhanced_counts[valid_genes, celltype_cells] <- enhanced_expr
    }
  }
  
  # 创建增强的SingleCellExperiment对象
  enhanced_sce <- SingleCellExperiment(
    assays = list(counts = enhanced_counts),
    colData = colData(sim),
    rowData = module_design$gene_assignments
  )
  
  # 添加logcounts和降维
  enhanced_sce <- logNormCounts(enhanced_sce)
  enhanced_sce <- runPCA(enhanced_sce)
  
  if (verbose) {
    log_message("增强单细胞数据生成完成:")
    log_message("  - 细胞数: %d", ncol(enhanced_sce))
    log_message("  - 基因数: %d", nrow(enhanced_sce))
    log_message("  - 细胞类型: %s", paste(levels(enhanced_sce$CellType), collapse = ", "))
  }
  
  return(enhanced_sce)
}

# 辅助函数: 应用共表达结构 --------------------------------------------------
#' @description 在基因表达矩阵中应用共表达结构
#' @param expr_matrix 表达矩阵 (基因 x 细胞)
#' @param correlation_strength 相关性强度 (0-1)
#' @param method 方法 ("factor_model" 或 "direct_correlation")
#' @return 调整后的表达矩阵
apply_coexpression_structure <- function(
    expr_matrix, 
    correlation_strength, 
    method = "factor_model"
) {
  if (nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    return(expr_matrix)
  }
  
  if (method == "factor_model") {
    # 使用因子模型方法，更符合生物学原理
    n_genes <- nrow(expr_matrix)
    n_cells <- ncol(expr_matrix)
    
    # 生成共同因子 (共同的调控信号)
    common_factor <- rnorm(n_cells, 0, correlation_strength)
    
    # 为每个基因生成因子载荷
    factor_loadings <- runif(n_genes, 0.3, 1.0) * 
      ifelse(runif(n_genes) > 0.2, 1, -1)  # 20%负相关
    
    # 应用因子结构
    for (i in 1:n_genes) {
      # 保留原始表达的一部分，添加共表达信号
      original_expr <- as.numeric(expr_matrix[i, ])
      factor_signal <- factor_loadings[i] * common_factor
      
      # 混合原始信号和共表达信号
      mixed_expr <- (1 - correlation_strength) * original_expr + 
        correlation_strength * factor_signal
      
      # 确保非负性
      expr_matrix[i, ] <- pmax(0, mixed_expr)
    }
    
  } else if (method == "direct_correlation") {
    # 直接相关性方法（较简单）
    # 计算基因间的目标相关性矩阵
    n_genes <- nrow(expr_matrix)
    target_corr <- matrix(correlation_strength, n_genes, n_genes)
    diag(target_corr) <- 1
    
    # 对表达矩阵进行Cholesky分解调整
    current_expr <- t(apply(expr_matrix, 1, scale))  # 标准化
    
    # 计算当前相关性
    current_corr <- cor(t(current_expr))
    current_corr[is.na(current_corr)] <- 0
    
    # 使用权重混合
    weight <- correlation_strength
    adjusted_corr <- weight * target_corr + (1 - weight) * current_corr
    
    # 确保矩阵为正定
    eigen_decomp <- eigen(adjusted_corr)
    adjusted_corr <- eigen_decomp$vectors %*% 
      diag(pmax(eigen_decomp$values, 0.01)) %*% 
      t(eigen_decomp$vectors)
    
    # 应用新的相关性结构
    chol_adj <- chol(adjusted_corr)
    new_expr <- chol_adj %*% current_expr
    
    # 重新缩放到原始范围
    for (i in 1:n_genes) {
      original_range <- range(expr_matrix[i, ])
      new_range <- range(new_expr[i, ])
      if (diff(new_range) > 0) {
        expr_matrix[i, ] <- (new_expr[i, ] - new_range[1]) / diff(new_range) * 
          diff(original_range) + original_range[1]
      }
    }
  }
  
  # 确保结果为非负整数
  expr_matrix[expr_matrix < 0] <- 0
  expr_matrix <- round(expr_matrix)
  
  return(expr_matrix)
}

# 模块3: 高级Pseudo-Bulk聚合器 (黄金标准实现) ------------------------------
#' @description 使用多种抽样策略和可变测序深度生成pseudo-bulk样本
#' @param sce 单细胞实验对象
#' @param n_samples 要生成的样本数量
#' @param sampling_strategy 抽样策略 ("random_mixture", "balanced", "unbalanced")
#' @param confounding_config 混杂配置（仅用于random_mixture）
#' @param sequencing_depth_config 测序深度配置
#' @param batch_config 批次配置
#' @param celltype_col 细胞类型列名
#' @param seed 随机种子
#' @param verbose 是否打印详细信息
#' @return 包含pseudo-bulk数据和详细元数据的列表
aggregate_to_pseudobulk_advanced <- function(
    sce,
    n_samples = 50,
    sampling_strategy = "random_mixture",
    confounding_config = list(
      level = "medium",  # "low", "medium", "high"
      proportion_ranges = NULL  # 自定义比例范围
    ),
    sequencing_depth_config = list(
      mean_depth = 500,     # 平均每样本细胞数
      depth_variation = 0.3, # 变异系数
      min_depth = 100,      # 最小细胞数
      max_depth = 1000      # 最大细胞数
    ),
    batch_config = list(
      n_batches = 3,
      batch_effect_strength = 0.5,
      confound_with_groups = TRUE  # 批次是否与生物学组别混杂
    ),
    celltype_col = "CellType",
    seed = 42,
    verbose = TRUE
) {
  set.seed(seed)
  
  if (verbose) {
    log_message("开始高级pseudo-bulk聚合...")
    log_message("抽样策略: %s", sampling_strategy)
    log_message("混杂水平: %s", confounding_config$level)
    log_message("平均测序深度: %d细胞", sequencing_depth_config$mean_depth)
  }
  
  # 获取细胞类型信息
  cell_types <- levels(sce[[celltype_col]])
  n_cell_types <- length(cell_types)
  
  # 为每个细胞类型创建细胞索引
  celltype_indices <- lapply(cell_types, function(ct) {
    which(sce[[celltype_col]] == ct)
  })
  names(celltype_indices) <- cell_types
  
  # 创建样本分组（两个主要生物学组）
  group_labels <- rep(c("GroupA", "GroupB"), length.out = n_samples)
  group_labels <- sample(group_labels)  # 随机化顺序
  
  # 创建批次分配
  batch_labels <- create_batch_assignment(
    n_samples = n_samples,
    group_labels = group_labels,
    batch_config = batch_config
  )
  
  # 初始化结果容器
  sample_metadata <- data.frame(
    sample_id = paste0("PseudoBulk_", sprintf("%03d", 1:n_samples)),
    biological_group = group_labels,
    batch = batch_labels,
    stringsAsFactors = FALSE
  )
  
  # 为每个样本生成细胞类型比例
  cell_proportions <- generate_cell_proportions(
    n_samples = n_samples,
    cell_types = cell_types,
    group_labels = group_labels,
    sampling_strategy = sampling_strategy,
    confounding_config = confounding_config,
    verbose = verbose
  )
  
  # 为每个样本确定测序深度（细胞数量）
  sequencing_depths <- generate_sequencing_depths(
    n_samples = n_samples,
    depth_config = sequencing_depth_config,
    verbose = verbose
  )
  
  # 初始化表达矩阵
  pseudobulk_counts <- matrix(0, nrow = nrow(sce), ncol = n_samples)
  rownames(pseudobulk_counts) <- rownames(sce)
  colnames(pseudobulk_counts) <- sample_metadata$sample_id
  
  # 存储实际的细胞组成信息
  actual_cell_counts <- matrix(0, nrow = n_samples, ncol = n_cell_types)
  colnames(actual_cell_counts) <- cell_types
  rownames(actual_cell_counts) <- sample_metadata$sample_id
  
  # 为每个样本进行细胞抽样和表达聚合
  if (verbose) log_message("开始为每个样本抽样细胞...")
  
  for (i in 1:n_samples) {
    if (verbose && i %% 10 == 0) {
      log_message("处理样本 %d/%d", i, n_samples)
    }
    
    sample_depth <- sequencing_depths[i]
    sample_proportions <- cell_proportions[i, ]
    
    # 根据比例和深度确定每种细胞类型的细胞数量
    target_cell_counts <- round(sample_depth * sample_proportions)
    
    # 确保至少有一些细胞
    target_cell_counts <- pmax(1, target_cell_counts)
    
    # 如果总数超过了目标深度，按比例缩减
    if (sum(target_cell_counts) > sample_depth) {
      scaling_factor <- sample_depth / sum(target_cell_counts)
      target_cell_counts <- pmax(1, round(target_cell_counts * scaling_factor))
    }
    
    # 抽样细胞并聚合表达
    selected_cells <- c()
    actual_counts <- numeric(n_cell_types)
    names(actual_counts) <- cell_types
    
    for (j in 1:n_cell_types) {
      cell_type <- cell_types[j]
      n_cells_needed <- target_cell_counts[j]
      available_cells <- celltype_indices[[cell_type]]
      
      if (length(available_cells) == 0) {
        actual_counts[j] <- 0
        next
      }
      
      # 抽样细胞（有放回或无放回取决于可用数量）
      if (length(available_cells) >= n_cells_needed) {
        sampled_cells <- sample(available_cells, n_cells_needed, replace = FALSE)
      } else {
        sampled_cells <- sample(available_cells, n_cells_needed, replace = TRUE)
      }
      
      selected_cells <- c(selected_cells, sampled_cells)
      actual_counts[j] <- length(sampled_cells)
    }
    
    # 聚合表达量
    if (length(selected_cells) > 0) {
      pseudobulk_counts[, i] <- rowSums(counts(sce)[, selected_cells, drop = FALSE])
    }
    
    # 存储实际细胞数量
    actual_cell_counts[i, ] <- actual_counts
  }
  
  # 应用批次效应
  if (verbose) log_message("应用批次效应...")
  pseudobulk_counts <- apply_batch_effects(
    count_matrix = pseudobulk_counts,
    batch_labels = batch_labels,
    batch_config = batch_config
  )
  
  # 计算质量控制指标
  sample_metadata$total_counts <- colSums(pseudobulk_counts)
  sample_metadata$detected_genes <- colSums(pseudobulk_counts > 0)
  sample_metadata$actual_cell_count <- rowSums(actual_cell_counts)
  
  # 添加细胞类型比例信息
  for (j in 1:n_cell_types) {
    # 预期比例
    sample_metadata[paste0("expected_prop_", cell_types[j])] <- cell_proportions[, j]
    # 实际比例
    sample_metadata[paste0("actual_prop_", cell_types[j])] <- 
      actual_cell_counts[, j] / rowSums(actual_cell_counts)
  }
  
  # 准备基因元数据
  gene_metadata <- as.data.frame(rowData(sce))
  if (nrow(gene_metadata) == 0) {
    gene_metadata <- data.frame(
      gene_id = rownames(sce),
      stringsAsFactors = FALSE
    )
  }
  
  # 汇总统计
  if (verbose) {
    log_message("Pseudo-bulk聚合完成:")
    log_message("  - 生成样本数: %d", n_samples)
    log_message("  - 基因数: %d", nrow(pseudobulk_counts))
    log_message("  - 平均文库大小: %.0f", mean(sample_metadata$total_counts))
    log_message("  - 文库大小范围: %.0f - %.0f", 
                min(sample_metadata$total_counts), 
                max(sample_metadata$total_counts))
    log_message("  - 平均检测基因数: %.0f", mean(sample_metadata$detected_genes))
  }
  
  # 返回结果
  return(list(
    counts = pseudobulk_counts,
    sample_metadata = sample_metadata,
    gene_metadata = gene_metadata,
    cell_counts = actual_cell_counts,
    design_info = list(
      sampling_strategy = sampling_strategy,
      confounding_config = confounding_config,
      sequencing_depth_config = sequencing_depth_config,
      batch_config = batch_config,
      cell_types = cell_types
    )
  ))
}

# 辅助函数: 批次分配 -------------------------------------------------------
#' @description 创建平衡或混杂的批次分配
create_batch_assignment <- function(n_samples, group_labels, batch_config) {
  n_batches <- batch_config$n_batches
  confound_with_groups <- batch_config$confound_with_groups
  
  if (confound_with_groups) {
    # 创建批次与组别的混杂
    unique_groups <- unique(group_labels)
    batch_assignment <- character(n_samples)
    
    # 为每个组优先分配到某些批次
    for (i in seq_along(unique_groups)) {
      group_samples <- which(group_labels == unique_groups[i])
      
      # 该组偏向的批次
      preferred_batches <- ((i - 1) %% n_batches) + 1
      
      # 80%分配到偏向批次，20%随机分配
      n_preferred <- round(0.8 * length(group_samples))
      preferred_indices <- sample(group_samples, n_preferred)
      remaining_indices <- setdiff(group_samples, preferred_indices)
      
      batch_assignment[preferred_indices] <- paste0("Batch", preferred_batches)
      if (length(remaining_indices) > 0) {
        batch_assignment[remaining_indices] <- 
          sample(paste0("Batch", 1:n_batches), length(remaining_indices), replace = TRUE)
      }
    }
  } else {
    # 平衡的批次分配
    batch_assignment <- rep(paste0("Batch", 1:n_batches), length.out = n_samples)
    batch_assignment <- sample(batch_assignment)
  }
  
  return(batch_assignment)
}

# 辅助函数: 细胞比例生成 ---------------------------------------------------
#' @description 根据不同策略生成细胞类型比例
generate_cell_proportions <- function(
    n_samples, 
    cell_types, 
    group_labels, 
    sampling_strategy, 
    confounding_config,
    verbose = TRUE
) {
  n_cell_types <- length(cell_types)
  
  if (sampling_strategy == "balanced") {
    # 平衡抽样：所有样本使用相同比例
    if (verbose) log_message("使用平衡抽样策略")
    base_proportions <- rep(1/n_cell_types, n_cell_types)
    proportions_matrix <- matrix(base_proportions, nrow = n_samples, ncol = n_cell_types, byrow = TRUE)
    
  } else if (sampling_strategy == "unbalanced") {
    # 非平衡抽样：固定数量策略
    if (verbose) log_message("使用非平衡抽样策略")
    # 为每种细胞类型分配固定的目标细胞数，然后归一化
    fixed_counts <- sample(50:200, n_cell_types, replace = TRUE)
    base_proportions <- fixed_counts / sum(fixed_counts)
    proportions_matrix <- matrix(base_proportions, nrow = n_samples, ncol = n_cell_types, byrow = TRUE)
    
  } else if (sampling_strategy == "random_mixture") {
    # 随机混合：可控混杂水平
    if (verbose) log_message("使用随机混合抽样策略，混杂水平: %s", confounding_config$level)
    
    # 定义混杂水平对应的比例变异范围
    variation_ranges <- switch(
      confounding_config$level,
      "low" = list(min_prop = 0.4, max_prop = 0.6),      # 低混杂：40%-60%
      "medium" = list(min_prop = 0.2, max_prop = 0.8),   # 中混杂：20%-80%
      "high" = list(min_prop = 0.05, max_prop = 0.95),   # 高混杂：5%-95%
      list(min_prop = 0.2, max_prop = 0.8)              # 默认中等
    )
    
    # 如果用户提供了自定义范围，使用自定义范围
    if (!is.null(confounding_config$proportion_ranges)) {
      variation_ranges <- confounding_config$proportion_ranges
    }
    
    proportions_matrix <- matrix(0, nrow = n_samples, ncol = n_cell_types)
    
    # 为每个样本生成随机比例
    for (i in 1:n_samples) {
      # 生成随机比例，确保在指定范围内
      sample_props <- runif(n_cell_types, 
                            variation_ranges$min_prop / n_cell_types,
                            variation_ranges$max_prop / n_cell_types)
      
      # 添加一些组别相关的偏向
      group <- group_labels[i]
      if (group == "GroupA") {
        # GroupA偏向前几种细胞类型
        bias_indices <- 1:ceiling(n_cell_types/2)
        sample_props[bias_indices] <- sample_props[bias_indices] * 1.5
      } else {
        # GroupB偏向后几种细胞类型
        bias_indices <- ceiling(n_cell_types/2 + 1):n_cell_types
        sample_props[bias_indices] <- sample_props[bias_indices] * 1.5
      }
      
      # 归一化
      sample_props <- sample_props / sum(sample_props)
      
      # 确保在合理范围内
      sample_props <- pmax(0.01, pmin(0.8, sample_props))
      sample_props <- sample_props / sum(sample_props)
      
      proportions_matrix[i, ] <- sample_props
    }
    
  } else {
    stop("不支持的抽样策略: ", sampling_strategy)
  }
  
  # 设置列名
  colnames(proportions_matrix) <- cell_types
  rownames(proportions_matrix) <- paste0("Sample_", 1:n_samples)
  
  return(proportions_matrix)
}

# 辅助函数: 测序深度生成 ---------------------------------------------------
#' @description 生成可变的测序深度（细胞数量）
generate_sequencing_depths <- function(n_samples, depth_config, verbose = TRUE) {
  mean_depth <- depth_config$mean_depth
  depth_variation <- depth_config$depth_variation
  min_depth <- depth_config$min_depth
  max_depth <- depth_config$max_depth
  
  if (verbose) {
    log_message("生成可变测序深度: 均值=%d, CV=%.2f", mean_depth, depth_variation)
  }
  
  # 使用gamma分布生成深度，确保正值
  shape <- 1 / (depth_variation^2)
  scale <- mean_depth / shape
  
  depths <- rgamma(n_samples, shape = shape, scale = scale)
  
  # 限制在指定范围内
  depths <- pmax(min_depth, pmin(max_depth, round(depths)))
  
  if (verbose) {
    log_message("实际深度范围: %d - %d (均值: %.1f)", 
                min(depths), max(depths), mean(depths))
  }
  
  return(depths)
}

# 辅助函数: 批次效应应用 ---------------------------------------------------
#' @description 应用真实的批次效应
apply_batch_effects <- function(count_matrix, batch_labels, batch_config) {
  n_genes <- nrow(count_matrix)
  n_samples <- ncol(count_matrix)
  batches <- unique(batch_labels)
  n_batches <- length(batches)
  effect_strength <- batch_config$batch_effect_strength
  
  # 选择批次敏感基因（60%的基因受批次影响）
  batch_sensitive_prop <- 0.6
  n_sensitive_genes <- round(n_genes * batch_sensitive_prop)
  sensitive_genes <- sample(1:n_genes, n_sensitive_genes)
  
  # 为每个批次创建基因特异性效应
  batch_effects <- array(1, dim = c(n_genes, n_batches))
  
  for (b in 1:n_batches) {
    # 为这个批次的敏感基因创建效应
    batch_effects[sensitive_genes, b] <- rlnorm(
      n_sensitive_genes, 
      meanlog = 0, 
      sdlog = effect_strength * 0.5
    )
  }
  
  # 应用批次效应
  adjusted_counts <- count_matrix
  for (i in 1:n_samples) {
    batch_idx <- which(batches == batch_labels[i])
    adjusted_counts[, i] <- adjusted_counts[, i] * batch_effects[, batch_idx]
  }
  
  # 确保结果为非负整数
  adjusted_counts[adjusted_counts < 0] <- 0
  adjusted_counts <- round(adjusted_counts)
  
  return(adjusted_counts)
}

# 模块4: 主控制函数 (参数化设计) -------------------------------------------
#' @description 运行完整的黄金标准pseudo-bulk生成管道
#' @param pipeline_config 管道配置参数
#' @param output_dir 输出目录
#' @param save_intermediate 是否保存中间结果
#' @param verbose 是否显示详细信息
#' @return 生成的pseudo-bulk数据和相关信息
run_golden_standard_pipeline <- function(
    pipeline_config = list(),
    output_dir = "golden_standard_results",
    save_intermediate = FALSE,
    verbose = TRUE
) {
  
  # 默认配置
  default_config <- list(
    # 模块设计参数
    module_design = list(
      n_cell_types = 6,
      n_genes_total = 3000,
      modules_per_celltype = 2,
      genes_per_module = c(30, 80),
      shared_modules_fraction = 0.2,
      correlation_strength = c(0.6, 0.9)
    ),
    
    # 单细胞生成参数
    single_cell = list(
      n_cells = 3000,
      group_probabilities = NULL,
      de_gene_proportion = 0.2
    ),
    
    # Pseudo-bulk聚合参数
    pseudobulk = list(
      n_samples = 50,
      sampling_strategy = "random_mixture",  # "random_mixture", "balanced", "unbalanced"
      confounding_config = list(
        level = "medium",  # "low", "medium", "high"
        proportion_ranges = NULL
      ),
      sequencing_depth_config = list(
        mean_depth = 500,
        depth_variation = 0.3,
        min_depth = 100,
        max_depth = 1000
      ),
      batch_config = list(
        n_batches = 3,
        batch_effect_strength = 0.5,
        confound_with_groups = TRUE
      )
    ),
    
    # 全局参数
    seed = 42
  )
  
  # 合并用户配置和默认配置
  config <- modifyList(default_config, pipeline_config)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (verbose) {
    log_message("开始黄金标准Pseudo-Bulk数据生成管道")
    log_message("配置摘要:")
    log_message("  - 细胞类型数: %d", config$module_design$n_cell_types)
    log_message("  - 基因数: %d", config$module_design$n_genes_total)
    log_message("  - 单细胞数: %d", config$single_cell$n_cells)
    log_message("  - Pseudo-bulk样本数: %d", config$pseudobulk$n_samples)
    log_message("  - 抽样策略: %s", config$pseudobulk$sampling_strategy)
  }
  
  # 第一步：设计细胞类型特异性共表达模块
  if (verbose) log_message("=== 第1步：设计共表达模块 ===")
  module_design <- do.call(design_celltype_specific_modules, 
                           c(config$module_design, seed = config$seed))
  
  # 第二步：生成增强的单细胞数据
  if (verbose) log_message("=== 第2步：生成增强单细胞数据 ===")
  enhanced_sce <- do.call(generate_enhanced_single_cell_data, 
                          c(list(module_design = module_design), 
                            config$single_cell, 
                            seed = config$seed,
                            verbose = verbose))
  
  # 第三步：聚合为pseudo-bulk样本
  if (verbose) log_message("=== 第3步：聚合为Pseudo-Bulk ===")
  pseudobulk_data <- do.call(aggregate_to_pseudobulk_advanced,
                             c(list(sce = enhanced_sce),
                               config$pseudobulk,
                               seed = config$seed,
                               verbose = verbose))
  
  # 第四步：保存结果
  if (verbose) log_message("=== 第4步：保存结果 ===")
  
  # 主要结果
  main_output_file <- file.path(output_dir, "pseudobulk_data.rds")
  saveRDS(pseudobulk_data, main_output_file)
  
  # 保存配置
  config_file <- file.path(output_dir, "pipeline_config.rds")
  saveRDS(config, config_file)
  
  # 可选：保存中间结果
  if (save_intermediate) {
    module_file <- file.path(output_dir, "module_design.rds")
    saveRDS(module_design, module_file)
    
    sce_file <- file.path(output_dir, "enhanced_single_cell.rds")
    saveRDS(enhanced_sce, sce_file)
  }
  
  # 生成质量控制报告
  if (verbose) log_message("=== 生成质量控制报告 ===")
  qc_report <- generate_qc_report(
    pseudobulk_data = pseudobulk_data,
    module_design = module_design,
    config = config,
    output_dir = output_dir
  )
  
  if (verbose) {
    log_message("管道执行完成!")
    log_message("主要输出文件: %s", main_output_file)
    log_message("质量控制报告: %s", file.path(output_dir, "qc_report.html"))
  }
  
  # 返回结果
  return(list(
    pseudobulk_data = pseudobulk_data,
    module_design = module_design,
    enhanced_sce = if(save_intermediate) enhanced_sce else NULL,
    config = config,
    output_files = list(
      main = main_output_file,
      config = config_file,
      qc_report = file.path(output_dir, "qc_report.html")
    )
  ))
}

# 模块5: 质量控制和验证 ----------------------------------------------------
#' @description 生成综合的质量控制报告
generate_qc_report <- function(pseudobulk_data, module_design, config, output_dir) {
  if (!require(rmarkdown, quietly = TRUE)) {
    warning("rmarkdown包未安装，跳过HTML报告生成")
    return(NULL)
  }
  
  # 创建简化的文本报告
  report_file <- file.path(output_dir, "qc_summary.txt")
  
  sink(report_file)
  cat("=== 黄金标准Pseudo-Bulk数据质量控制报告 ===\n")
  cat("生成时间:", format(Sys.time()), "\n\n")
  
  # 基本统计
  cat("## 基本统计信息\n")
  cat("样本数量:", ncol(pseudobulk_data$counts), "\n")
  cat("基因数量:", nrow(pseudobulk_data$counts), "\n")
  cat("细胞类型数:", length(config$module_design$n_cell_types), "\n")
  
  # 文库大小统计
  cat("\n## 文库大小分布\n")
  lib_sizes <- pseudobulk_data$sample_metadata$total_counts
  cat("平均文库大小:", round(mean(lib_sizes)), "\n")
  cat("文库大小范围:", min(lib_sizes), "-", max(lib_sizes), "\n")
  cat("文库大小CV:", round(sd(lib_sizes)/mean(lib_sizes), 3), "\n")
  
  # 细胞类型比例统计
  cat("\n## 细胞类型比例\n")
  prop_cols <- grep("^actual_prop_", colnames(pseudobulk_data$sample_metadata), value = TRUE)
  for (col in prop_cols) {
    celltype <- gsub("^actual_prop_", "", col)
    props <- pseudobulk_data$sample_metadata[[col]]
    cat(sprintf("%s: %.3f ± %.3f (%.3f-%.3f)\n", 
                celltype, mean(props), sd(props), min(props), max(props)))
  }
  
  # 批次分布
  cat("\n## 批次分布\n")
  batch_table <- table(pseudobulk_data$sample_metadata$batch, 
                       pseudobulk_data$sample_metadata$biological_group)
  print(batch_table)
  
  cat("\n## 共表达模块统计\n")
  for (celltype in names(module_design$celltype_modules)) {
    n_modules <- length(module_design$celltype_modules[[celltype]])
    cat(sprintf("%s: %d个模块\n", celltype, n_modules))
  }
  
  sink()
  
  log_message("质量控制摘要已保存至: %s", report_file)
  
  return(report_file)
}

# 模块6: 便捷接口函数 ------------------------------------------------------
#' @description 便捷的主函数，提供预设配置选项
#' @param preset 预设配置 ("simple", "complex", "high_confounding", "custom")
#' @param custom_config 自定义配置（当preset="custom"时使用）
#' @param output_dir 输出目录
#' @param verbose 是否显示详细信息
#' @return 生成的数据和报告
generate_pseudobulk_data <- function(
    preset = "simple",
    custom_config = NULL,
    output_dir = "pseudobulk_output",
    verbose = TRUE
) {
  
  # 预设配置
  preset_configs <- list(
    simple = list(
      module_design = list(
        n_cell_types = 4,
        n_genes_total = 2000,
        modules_per_celltype = 1,
        genes_per_module = c(20, 50)
      ),
      single_cell = list(n_cells = 1500),
      pseudobulk = list(
        n_samples = 30,
        sampling_strategy = "random_mixture",
        confounding_config = list(level = "medium")
      )
    ),
    
    complex = list(
      module_design = list(
        n_cell_types = 8,
        n_genes_total = 5000,
        modules_per_celltype = 3,
        genes_per_module = c(40, 100),
        shared_modules_fraction = 0.3
      ),
      single_cell = list(n_cells = 5000),
      pseudobulk = list(
        n_samples = 80,
        sampling_strategy = "random_mixture",
        confounding_config = list(level = "medium")
      )
    ),
    
    high_confounding = list(
      module_design = list(
        n_cell_types = 6,
        n_genes_total = 3000
      ),
      pseudobulk = list(
        n_samples = 60,
        sampling_strategy = "random_mixture",
        confounding_config = list(level = "high"),
        batch_config = list(
          n_batches = 4,
          batch_effect_strength = 0.8,
          confound_with_groups = TRUE
        )
      )
    )
  )
  
  # 选择配置
  if (preset == "custom") {
    if (is.null(custom_config)) {
      stop("使用custom预设时必须提供custom_config参数")
    }
    final_config <- custom_config
  } else if (preset %in% names(preset_configs)) {
    final_config <- preset_configs[[preset]]
  } else {
    stop("不支持的预设: ", preset, ". 支持的预设: ", 
         paste(names(preset_configs), collapse = ", "), ", custom")
  }
  
  if (verbose) {
    log_message("使用预设配置: %s", preset)
  }
  
  # 运行管道
  result <- run_golden_standard_pipeline(
    pipeline_config = final_config,
    output_dir = output_dir,
    save_intermediate = FALSE,
    verbose = verbose
  )
  
  return(result)
}

# 示例使用函数 -------------------------------------------------------------
#' @description 提供使用示例
demo_usage <- function() {
  cat("=== 黄金标准Pseudo-Bulk数据生成器使用示例 ===\n\n")
  
  cat("1. 简单示例：\n")
  cat("result <- generate_pseudobulk_data(preset = 'simple', output_dir = 'demo_simple')\n\n")
  
  cat("2. 复杂场景：\n")
  cat("result <- generate_pseudobulk_data(preset = 'complex', output_dir = 'demo_complex')\n\n")
  
  cat("3. 高混杂场景：\n")
  cat("result <- generate_pseudobulk_data(preset = 'high_confounding', output_dir = 'demo_confound')\n\n")
  
  cat("4. 自定义配置：\n")
  cat("custom_cfg <- list(\n")
  cat("  module_design = list(n_cell_types = 10, n_genes_total = 4000),\n")
  cat("  pseudobulk = list(\n")
  cat("    sampling_strategy = 'balanced',\n")
  cat("    n_samples = 100,\n")
  cat("    sequencing_depth_config = list(mean_depth = 800, depth_variation = 0.4)\n")
  cat("  )\n")
  cat(")\n")
  cat("result <- generate_pseudobulk_data(preset = 'custom', custom_config = custom_cfg)\n\n")
  
  cat("生成的数据包含：\n")
  cat("- result$pseudobulk_data$counts: 表达矩阵\n")
  cat("- result$pseudobulk_data$sample_metadata: 样本元数据\n")
  cat("- result$pseudobulk_data$gene_metadata: 基因元数据\n")
  cat("- result$module_design: 共表达模块设计\n")
}

# 在脚本末尾显示使用示例
if (interactive()) {
  demo_usage()
}
