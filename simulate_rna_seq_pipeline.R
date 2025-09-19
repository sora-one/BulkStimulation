# 包加载
ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("正在安装包:", pkg))
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
}

# 加载必要的库
required_packages <- c(
  "BiocParallel", "splatter", "SingleCellExperiment", "scater", "Matrix", 
  "ggplot2", "dplyr", "gtools", "rmarkdown"
)
suppressPackageStartupMessages(ensure_packages(required_packages))


# --- 修复: 修正log_message函数签名以正确处理格式化参数 ---
log_message <- function(msg, ..., level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  dots <- list(...)
  if (length(dots) > 0) {
    tryCatch({
      formatted_msg <- sprintf(msg, ...)
    }, error = function(e) {
      warning("sprintf formatting failed for message '", msg, "': ", e$message)
      formatted_msg <<- paste(msg, "| Args:", paste(dots, collapse = ", "))
    })
  } else {
    formatted_msg <- msg
  }
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, formatted_msg))
}


# 模块1: 细胞类型特异性共表达模块设计器 (无变动)
design_celltype_specific_modules <- function(
    n_cell_types = 6, n_genes_total = 3000, modules_per_celltype = 2,
    genes_per_module = c(30, 80), shared_modules_fraction = 0.2,
    correlation_strength = c(0.6, 0.9), seed = 42
) {
  set.seed(seed)
  log_message("设计细胞类型特异性共表达模块...")
  gene_names <- paste0("Gene", sprintf("%04d", 1:n_genes_total))
  celltype_modules <- setNames(vector("list", n_cell_types), paste0("CellType", 1:n_cell_types))
  shared_modules <- list()
  gene_assignments <- data.frame(gene_id = gene_names, stringsAsFactors = FALSE)
  for (i in 1:n_cell_types) {
    gene_assignments[[paste0("in_CellType", i)]] <- FALSE
  }
  
  total_modules_target <- n_cell_types * modules_per_celltype
  n_shared_modules <- round(total_modules_target * shared_modules_fraction)
  n_specific_modules <- total_modules_target - n_shared_modules
  used_genes <- c()
  
  if (n_shared_modules > 0) {
    log_message("创建%d个跨细胞类型共享模块...", n_shared_modules)
    for (i in 1:n_shared_modules) {
      n_sharing_types <- sample(2:min(4, n_cell_types), 1)
      sharing_types <- sample(1:n_cell_types, n_sharing_types)
      module_size <- sample(genes_per_module[1]:genes_per_module[2], 1)
      available_genes <- setdiff(gene_names, used_genes)
      if (length(available_genes) < module_size) { if (length(available_genes) == 0) break; module_size <- length(available_genes) }
      module_genes <- sample(available_genes, module_size)
      used_genes <- c(used_genes, module_genes)
      corr_strength <- runif(1, correlation_strength[1], correlation_strength[2])
      module_name <- paste0("SharedModule", i)
      shared_modules[[module_name]] <- list(genes = module_genes, cell_types = paste0("CellType", sharing_types), correlation_strength = corr_strength, module_type = "shared")
      for (ct_idx in sharing_types) {
        celltype_name <- paste0("CellType", ct_idx)
        gene_assignments[gene_assignments$gene_id %in% module_genes, paste0("in_", celltype_name)] <- TRUE
        celltype_modules[[celltype_name]][[module_name]] <- list(genes = module_genes, correlation_strength = corr_strength, module_type = "shared")
      }
    }
  }
  
  log_message("创建细胞类型特异性模块...")
  specific_modules_per_type_base <- n_specific_modules %/% n_cell_types
  remaining_specific_modules <- n_specific_modules %% n_cell_types
  for (i in 1:n_cell_types) {
    celltype_name <- paste0("CellType", i)
    n_modules_this_type <- specific_modules_per_type_base + (i <= remaining_specific_modules)
    if (n_modules_this_type == 0) next
    for (j in 1:n_modules_this_type) {
      module_size <- sample(genes_per_module[1]:genes_per_module[2], 1)
      available_genes <- setdiff(gene_names, used_genes)
      if (length(available_genes) < module_size) { if (length(available_genes) == 0) break; module_size <- length(available_genes) }
      module_genes <- sample(available_genes, module_size)
      used_genes <- c(used_genes, module_genes)
      corr_strength <- runif(1, correlation_strength[1] + 0.1, correlation_strength[2])
      module_name <- paste0(celltype_name, "_SpecificModule", j)
      celltype_modules[[celltype_name]][[module_name]] <- list(genes = module_genes, correlation_strength = corr_strength, module_type = "specific")
      gene_assignments[gene_assignments$gene_id %in% module_genes, paste0("in_", celltype_name)] <- TRUE
    }
  }
  
  total_assigned_genes_count <- length(unique(used_genes))
  log_message("  - 总模块数: %d (共享: %d, 特异: %d)", n_shared_modules + n_specific_modules, n_shared_modules, n_specific_modules)
  log_message("  - 已分配基因: %d/%d (%.1f%%)", total_assigned_genes_count, n_genes_total, 100 * total_assigned_genes_count / n_genes_total)
  
  return(list(celltype_modules = celltype_modules, shared_modules = shared_modules, gene_assignments = gene_assignments, design_params = list(n_cell_types = n_cell_types, n_genes_total = n_genes_total, modules_per_celltype = modules_per_celltype, genes_per_module = genes_per_module, shared_modules_fraction = shared_modules_fraction, correlation_strength = correlation_strength)))
}

# 模块2: 增强的单细胞数据生成器 (无变动)
generate_enhanced_single_cell_data <- function(
    module_design, n_cells = 3000, group_probabilities = NULL,
    de_gene_proportion = 0.2, seed = 42, verbose = TRUE
) {
  set.seed(seed)
  if (verbose) log_message("生成具有预定义共表达结构的单细胞数据...")
  n_cell_types <- module_design$design_params$n_cell_types
  n_genes <- module_design$design_params$n_genes_total
  if (is.null(group_probabilities)) { group_probabilities <- rep(1/n_cell_types, n_cell_types) }
  params <- newSplatParams(batchCells = n_cells, nGenes = n_genes, group.prob = group_probabilities, de.prob = rep(de_gene_proportion, n_cell_types), de.facLoc = 1.0, seed = seed)
  if (verbose) log_message("使用Splatter生成基础单细胞表达谱...")
  sim <- splatSimulate(params, method = "groups", verbose = FALSE)
  cell_type_levels <- paste0("CellType", 1:n_cell_types)
  colData(sim)$CellType <- factor(gsub("Group", "CellType", sim$Group), levels = cell_type_levels)
  
  if (verbose) log_message("应用细胞类型特异性共表达结构...")
  enhanced_counts <- as(counts(sim), "sparseMatrix")
  for (celltype_name in names(module_design$celltype_modules)) {
    if (verbose) log_message("处理%s的共表达模块...", celltype_name)
    celltype_cells <- which(sim$CellType == celltype_name)
    if (length(celltype_cells) == 0) next
    celltype_modules <- module_design$celltype_modules[[celltype_name]]
    for (module_name in names(celltype_modules)) {
      module_info <- celltype_modules[[module_name]]
      valid_genes <- intersect(module_info$genes, rownames(enhanced_counts))
      if (length(valid_genes) < 2) next
      current_expr <- enhanced_counts[valid_genes, celltype_cells, drop = FALSE]
      enhanced_expr <- apply_coexpression_structure(expr_matrix = current_expr, correlation_strength = module_info$correlation_strength)
      enhanced_counts[valid_genes, celltype_cells] <- enhanced_expr
    }
  }
  
  enhanced_sce <- SingleCellExperiment(assays = list(counts = enhanced_counts), colData = colData(sim), rowData = module_design$gene_assignments)
  enhanced_sce <- logNormCounts(enhanced_sce)
  enhanced_sce <- runPCA(enhanced_sce)
  
  if (verbose) { log_message("增强单细胞数据生成完成: %d 细胞, %d 基因", ncol(enhanced_sce), nrow(enhanced_sce)) }
  return(enhanced_sce)
}

# 辅助函数: 应用共表达结构 (无变动)
apply_coexpression_structure <- function(expr_matrix, correlation_strength) {
  if (nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) return(expr_matrix)
  n_genes <- nrow(expr_matrix)
  n_cells <- ncol(expr_matrix)
  log_expr <- log1p(as.matrix(expr_matrix))
  common_factor <- rnorm(n_cells, mean = 0, sd = 1)
  factor_loadings <- runif(n_genes, 0.5, 1.5) * sample(c(1, -1), n_genes, replace = TRUE, prob = c(0.8, 0.2))
  mixed_log_expr <- matrix(0, nrow = n_genes, ncol = n_cells)
  for (i in 1:n_genes) {
    mixed_log_expr[i, ] <- sqrt(1 - correlation_strength^2) * log_expr[i, ] + 
      correlation_strength * factor_loadings[i] * common_factor
  }
  enhanced_expr <- expm1(mixed_log_expr)
  enhanced_expr[enhanced_expr < 0] <- 0
  return(round(enhanced_expr))
}

# --- 增强: aggregate_to_pseudobulk_advanced 现在接受 'cell_proportions' 参数 ---
aggregate_to_pseudobulk_advanced <- function(
    sce, n_samples = 50, sampling_strategy = "random_mixture",
    cell_proportions = NULL, # 新增参数！
    confounding_config = list(level = "medium"),
    sequencing_depth_config = list(mean_depth = 500, depth_variation = 0.3, min_depth = 100, max_depth = 1000),
    batch_config = list(n_batches = 3, batch_effect_strength = 0.5, confound_with_groups = TRUE),
    celltype_col = "CellType", seed = 42, verbose = TRUE
) {
  set.seed(seed)
  if (verbose) log_message("开始高级pseudo-bulk聚合...")
  
  cell_types <- levels(sce[[celltype_col]])
  n_cell_types <- length(cell_types)
  celltype_indices <- lapply(cell_types, function(ct) which(sce[[celltype_col]] == ct))
  names(celltype_indices) <- cell_types
  
  group_labels <- sample(rep(c("GroupA", "GroupB"), length.out = n_samples))
  batch_labels <- create_batch_assignment(n_samples, group_labels, batch_config)
  
  sample_metadata <- data.frame(sample_id = paste0("PseudoBulk_", sprintf("%03d", 1:n_samples)), biological_group = group_labels, batch = batch_labels, stringsAsFactors = FALSE)
  
  # --- 增强逻辑: 如果提供了cell_proportions，则使用它；否则，根据策略生成 ---
  if (is.null(cell_proportions)) {
    log_message("根据策略 '%s' 生成细胞比例...", sampling_strategy)
    cell_proportions <- generate_cell_proportions(n_samples, cell_types, group_labels, sampling_strategy, confounding_config, verbose)
  } else {
    log_message("使用用户提供的自定义细胞比例矩阵。")
    # 确保提供的矩阵维度正确
    stopifnot(nrow(cell_proportions) == n_samples)
    stopifnot(ncol(cell_proportions) == n_cell_types)
  }
  
  sequencing_depths <- generate_sequencing_depths(n_samples, sequencing_depth_config, verbose)
  
  pseudobulk_counts <- Matrix::Matrix(0, nrow = nrow(sce), ncol = n_samples, sparse = TRUE, dimnames = list(rownames(sce), sample_metadata$sample_id))
  actual_cell_counts <- matrix(0, nrow = n_samples, ncol = n_cell_types, dimnames = list(sample_metadata$sample_id, cell_types))
  
  sc_counts_matrix <- counts(sce)
  
  if (verbose) log_message("为%d个样本抽样细胞...", n_samples)
  for (i in 1:n_samples) {
    sample_depth <- sequencing_depths[i]
    target_cell_counts <- round(sample_depth * cell_proportions[i, ])
    current_sum <- sum(target_cell_counts)
    if (current_sum > 0) { target_cell_counts <- round(target_cell_counts * (sample_depth / current_sum)) }
    if (sum(target_cell_counts) == 0 && sample_depth > 0) { target_cell_counts[sample.int(n_cell_types, 1)] <- 1 }
    
    selected_cells <- unlist(lapply(seq_along(cell_types), function(j) {
      n_cells_needed <- target_cell_counts[j]
      if (n_cells_needed == 0) return(NULL)
      available_cells <- celltype_indices[[j]]
      if (length(available_cells) == 0) return(NULL)
      replace_sampling <- length(available_cells) < n_cells_needed
      sampled <- sample(available_cells, n_cells_needed, replace = replace_sampling)
      actual_cell_counts[i, j] <<- length(sampled)
      return(sampled)
    }))
    
    if (length(selected_cells) > 0) {
      pseudobulk_counts[, i] <- Matrix::rowSums(sc_counts_matrix[, selected_cells, drop = FALSE])
    }
  }
  
  if (verbose) log_message("应用批次效应...")
  pseudobulk_counts <- apply_batch_effects(pseudobulk_counts, batch_labels, batch_config)
  
  sample_metadata$total_counts <- Matrix::colSums(pseudobulk_counts)
  sample_metadata$detected_genes <- Matrix::colSums(pseudobulk_counts > 0)
  sample_metadata$actual_cell_count <- rowSums(actual_cell_counts)
  
  actual_proportions <- actual_cell_counts / pmax(1, rowSums(actual_cell_counts))
  colnames(cell_proportions) <- paste0("expected_prop_", colnames(cell_proportions))
  colnames(actual_proportions) <- paste0("actual_prop_", colnames(actual_proportions))
  sample_metadata <- cbind(sample_metadata, as.data.frame(cell_proportions), as.data.frame(actual_proportions))
  
  if (verbose) log_message("Pseudo-bulk聚合完成。")
  return(list(counts = pseudobulk_counts, sample_metadata = sample_metadata, gene_metadata = as.data.frame(rowData(sce)), cell_counts = actual_cell_counts, design_info = list(sampling_strategy = sampling_strategy, confounding_config = confounding_config, sequencing_depth_config = sequencing_depth_config, batch_config = batch_config, cell_types = cell_types)))
}


# 辅助函数: 批次分配，细胞比例生成等 (无变动)
create_batch_assignment <- function(n_samples, group_labels, batch_config) {
  n_batches <- batch_config$n_batches
  if (batch_config$confound_with_groups) {
    unique_groups <- unique(group_labels)
    batch_assignment <- character(n_samples)
    for (i in seq_along(unique_groups)) {
      group_samples <- which(group_labels == unique_groups[i])
      preferred_batch <- paste0("Batch", ((i - 1) %% n_batches) + 1)
      n_preferred <- min(length(group_samples), round(0.8 * length(group_samples)))
      preferred_indices <- sample(group_samples, n_preferred)
      batch_assignment[preferred_indices] <- preferred_batch
      remaining_indices <- setdiff(group_samples, preferred_indices)
      if (length(remaining_indices) > 0) { batch_assignment[remaining_indices] <- sample(paste0("Batch", 1:n_batches), length(remaining_indices), replace = TRUE) }
    }
  } else {
    batch_assignment <- sample(rep(paste0("Batch", 1:n_batches), length.out = n_samples))
  }
  return(batch_assignment)
}

generate_cell_proportions <- function(n_samples, cell_types, group_labels, sampling_strategy, confounding_config, verbose) {
  n_cell_types <- length(cell_types)
  proportions_matrix <- matrix(0, nrow = n_samples, ncol = n_cell_types)
  
  if (sampling_strategy == "balanced") {
    for(i in 1:n_samples) proportions_matrix[i,] <- gtools::rdirichlet(1, rep(10, n_cell_types))
  } else { # For "unbalanced" and "random_mixture"
    alpha_factor <- switch(confounding_config$level, "low" = 10, "medium" = 3, "high" = 0.5, 3)
    for (i in 1:n_samples) {
      alpha_base <- rep(1, n_cell_types)
      if (group_labels[i] == "GroupA") { alpha_base[1:ceiling(n_cell_types/2)] <- 3 } 
      else { alpha_base[(ceiling(n_cell_types/2) + 1):n_cell_types] <- 3 }
      proportions_matrix[i, ] <- gtools::rdirichlet(1, alpha_base * alpha_factor)
    }
  }
  colnames(proportions_matrix) <- cell_types
  return(proportions_matrix)
}

generate_sequencing_depths <- function(n_samples, depth_config, verbose) {
  if (depth_config$depth_variation <= 0) return(rep(depth_config$mean_depth, n_samples))
  shape <- 1 / (depth_config$depth_variation^2)
  scale <- depth_config$mean_depth / shape
  depths <- rgamma(n_samples, shape = shape, scale = scale)
  return(pmax(depth_config$min_depth, pmin(depth_config$max_depth, round(depths))))
}

apply_batch_effects <- function(count_matrix, batch_labels, batch_config) {
  n_genes <- nrow(count_matrix)
  batches <- unique(batch_labels)
  n_batches <- length(batches)
  sensitive_genes <- sample.int(n_genes, round(n_genes * 0.6))
  batch_multipliers <- matrix(1.0, nrow = n_genes, ncol = n_batches)
  for (b in 1:n_batches) {
    batch_multipliers[sensitive_genes, b] <- rlnorm(length(sensitive_genes), 0, batch_config$batch_effect_strength * 0.5)
  }
  batch_map <- match(batch_labels, batches)
  adjusted_counts <- count_matrix
  for(g in 1:n_genes){
    adjusted_counts[g,] <- adjusted_counts[g,] * batch_multipliers[g, batch_map]
  }
  adjusted_counts@x <- round(adjusted_counts@x)
  return(adjusted_counts)
}

# 模块4: 主控制函数 (无变动)
run_golden_standard_pipeline <- function(
    pipeline_config = list(), output_dir = "golden_standard_results",
    save_intermediate = FALSE, verbose = TRUE
) {
  default_config <- list(
    module_design = list(n_cell_types = 6, n_genes_total = 3000, modules_per_celltype = 2, genes_per_module = c(30, 80), shared_modules_fraction = 0.2, correlation_strength = c(0.6, 0.9)),
    single_cell = list(n_cells = 3000, de_gene_proportion = 0.2),
    pseudobulk = list(n_samples = 50, sampling_strategy = "random_mixture", confounding_config = list(level = "medium"), sequencing_depth_config = list(mean_depth = 500, depth_variation = 0.3, min_depth = 100, max_depth = 1000), batch_config = list(n_batches = 3, batch_effect_strength = 0.5, confound_with_groups = TRUE)),
    seed = 42
  )
  config <- modifyList(default_config, pipeline_config)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (verbose) log_message("开始黄金标准管道...")
  
  module_design <- do.call(design_celltype_specific_modules, c(config$module_design, seed = config$seed))
  enhanced_sce <- do.call(generate_enhanced_single_cell_data, c(list(module_design = module_design), config$single_cell, seed = config$seed, verbose = verbose))
  pseudobulk_data <- do.call(aggregate_to_pseudobulk_advanced, c(list(sce = enhanced_sce), config$pseudobulk, seed = config$seed, verbose = verbose))
  
  if (verbose) log_message("保存结果到 %s...", output_dir)
  saveRDS(pseudobulk_data, file.path(output_dir, "pseudobulk_data.rds"))
  saveRDS(config, file.path(output_dir, "pipeline_config.rds"))
  generate_qc_report(pseudobulk_data, module_design, config, output_dir)
  
  if (verbose) log_message("管道执行完成!")
  return(list(pseudobulk_data = pseudobulk_data, module_design = module_design, config = config))
}

# 模块5: 质量控制和验证 (无变动)
generate_qc_report <- function(pseudobulk_data, module_design, config, output_dir) {
  report_file <- file.path(output_dir, "qc_summary.txt")
  tryCatch({
    sink(report_file)
    cat("=== 黄金标准Pseudo-Bulk数据质量控制报告 ===\n\n")
    cat("## 基本统计\n")
    cat("样本数:", ncol(pseudobulk_data$counts), "\n")
    cat("基因数:", nrow(pseudobulk_data$counts), "\n")
    cat("细胞类型数:", config$module_design$n_cell_types, "\n\n")
    cat("## 文库大小分布\n")
    lib_sizes <- pseudobulk_data$sample_metadata$total_counts
    cat("平均文库大小:", round(mean(lib_sizes)), "\n")
    cat("文库大小范围:", min(lib_sizes), "-", max(lib_sizes), "\n\n")
    cat("## 细胞类型比例 (实际抽样)\n")
    prop_cols <- grep("^actual_prop_", colnames(pseudobulk_data$sample_metadata), value = TRUE)
    for (col in prop_cols) {
      props <- pseudobulk_data$sample_metadata[[col]]
      cat(sprintf("%s: %.3f ± %.3f\n", gsub("^actual_prop_", "", col), mean(props), sd(props)))
    }
    cat("\n## 批次与组别分布\n")
    print(table(Batch = pseudobulk_data$sample_metadata$batch, Group = pseudobulk_data$sample_metadata$biological_group))
    sink()
    log_message("质量控制摘要已保存至: %s", report_file)
  }, error = function(e) {
    if(sink.number() > 0) sink()
    log_message("生成QC报告失败: %s", e$message, level = "ERROR")
  })
  return(report_file)
}

# 模块6: 便捷接口函数 (无变动)
generate_pseudobulk_data <- function(
    preset = "simple", custom_config = NULL, output_dir = "pseudobulk_output", verbose = TRUE
) {
  preset_configs <- list(
    simple = list(
      module_design = list(n_cell_types = 4, n_genes_total = 2000, modules_per_celltype = 1, genes_per_module = c(20, 50)),
      single_cell = list(n_cells = 1500),
      pseudobulk = list(n_samples = 30, sampling_strategy = "random_mixture", confounding_config = list(level = "medium"), sequencing_depth_config = list(mean_depth = 300, depth_variation = 0.2, min_depth = 50, max_depth = 800), batch_config = list(n_batches = 2, batch_effect_strength = 0.3, confound_with_groups = FALSE))
    ),
    complex = list(
      module_design = list(n_cell_types = 8, n_genes_total = 5000, modules_per_celltype = 3, genes_per_module = c(40, 100), shared_modules_fraction = 0.3),
      single_cell = list(n_cells = 5000),
      pseudobulk = list(n_samples = 80, sampling_strategy = "random_mixture", confounding_config = list(level = "medium"))
    ),
    high_confounding = list(
      module_design = list(n_cell_types = 6, n_genes_total = 3000),
      pseudobulk = list(n_samples = 60, sampling_strategy = "random_mixture", confounding_config = list(level = "high"), batch_config = list(n_batches = 4, batch_effect_strength = 0.8, confound_with_groups = TRUE))
    ),
    balanced_composition = list(
      module_design = list(n_cell_types = 5, n_genes_total = 2500),
      single_cell = list(n_cells = 2000),
      pseudobulk = list(n_samples = 40, sampling_strategy = "balanced", confounding_config = list(level = "low"))
    )
  )
  
  if (preset == "custom") {
    if (is.null(custom_config)) stop("使用'custom'预设时必须提供'custom_config'参数")
    final_config <- custom_config
  } else if (preset %in% names(preset_configs)) {
    final_config <- preset_configs[[preset]]
  } else { stop("不支持的预设: ", preset) }
  
  if (verbose) log_message("使用预设配置: '%s'", preset)
  return(run_golden_standard_pipeline(pipeline_config = final_config, output_dir = output_dir, verbose = verbose))
}
# 示例使用函数 -------------------------------------------------------------
demo_usage <- function() {
  cat("=== 黄金标准Pseudo-Bulk数据生成器使用示例 ===\n\n")
  cat("1. 简单示例: > result <- generate_pseudobulk_data(preset = 'simple')\n")
  cat("2. 复杂场景: > result <- generate_pseudobulk_data(preset = 'complex')\n")
  cat("3. 高混杂场景: > result <- generate_pseudobulk_data(preset = 'high_confounding')\n")
  cat("4. 平衡组分: > result <- generate_pseudobulk_data(preset = 'balanced_composition')\n")
}

# 交互式会话中运行演示
if (interactive()) {
  demo_usage()
  cat("\n--- 正在运行一个简单的演示，请稍候... ---\n")
  # 确保gtools已安装
  ensure_packages("gtools")
  result_demo <- generate_pseudobulk_data(preset = "simple", output_dir = "simple_demo_output")
  cat("\n--- 演示完成！请查看 'simple_demo_output' 目录。---\n")
  print(head(result_demo$pseudobulk_data$sample_metadata))
}
