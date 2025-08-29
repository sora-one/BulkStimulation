# ============================================================================
# Confounding Factor Correction for RNA-Seq Data
# ============================================================================
# This script contains functions to apply various correction methods
# for known and unknown confounding factors in RNA-seq data, as described in
# Cote et al. (2022, DOI: 10.1186/s13059-022-02606-0)
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)       # For variance stabilizing transformation
  library(limma)        # For removeBatchEffect and voom
  library(sva)          # For SVA and num.sv functions
  library(matrixStats)  # For rowVars
  library(PEER)         # For PEER correction
  library(RUVcorr)      # For RUVCorr
  library(BiocParallel) # For parallel processing
  library(ggplot2)      # For visualization
  library(dplyr)        # For data manipulation
})

#' Apply Multiple Confounding Factor Correction Methods to RNA-Seq Data
#'
#' This function applies various correction methods to RNA-seq data to remove
#' the effects of technical and biological confounding factors.
#'
#' @param sim_data A list containing simulated RNA-seq data with elements:
#'   - bulk_counts: A matrix of raw gene expression counts (genes x samples)
#'   - bulk_metadata: A data frame with sample metadata
#'   - gene_metadata: A data frame with gene metadata
#'
#' @return A list of corrected expression matrices, each named after the
#'         correction method applied.
#'
apply_all_corrections <- function(sim_data) {
  
  # Extract data from the input list
  bulk_counts <- sim_data$bulk_counts
  bulk_metadata <- sim_data$bulk_metadata
  gene_metadata <- sim_data$gene_metadata
  
  # Check input data
  if (is.null(bulk_counts) || is.null(bulk_metadata) || is.null(gene_metadata)) {
    stop("Input data must contain bulk_counts, bulk_metadata, and gene_metadata")
  }
  
  # Make sure sample names match between counts and metadata
  if (!identical(colnames(bulk_counts), rownames(bulk_metadata))) {
    stop("Sample names in bulk_counts and bulk_metadata do not match")
  }
  
  # Make sure gene names match between counts and metadata
  if (!identical(rownames(bulk_counts), rownames(gene_metadata))) {
    stop("Gene names in bulk_counts and gene_metadata do not match")
  }
  
  # Create a list to store all corrected matrices
  corrected_matrices <- list()
  
  # ============================================================================
  # 1. Data Normalization
  # ============================================================================
  message("Performing variance stabilizing transformation...")
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = bulk_counts,
    colData = bulk_metadata,
    design = ~ group
  )
  
  # Apply variance stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  
  # Extract the normalized matrix
  normalized_matrix <- assay(vsd)
  
  # Add to the results list
  corrected_matrices$normalized_uncorrected <- normalized_matrix
  
  # ============================================================================
  # 2. Known Covariate Regression
  # ============================================================================
  message("Applying known covariate regression...")
  
  # Create design matrix for the biological factor (group)
  design <- model.matrix(~ group, data = bulk_metadata)
  
  # Identify cell type proportion columns (starting with "prop_")
  prop_cols <- grep("^prop_", names(bulk_metadata), value = TRUE)
  
  if (length(prop_cols) == 0) {
    warning("No cell type proportion columns found in bulk_metadata")
    covariates <- model.matrix(~ batch, data = bulk_metadata)[, -1, drop = FALSE]
  } else {
    # Create covariate matrix for batch and cell type proportions
    covariates <- model.matrix(
      as.formula(paste("~ batch +", paste(prop_cols, collapse = " + "))),
      data = bulk_metadata
    )[, -1, drop = FALSE]  # Remove intercept
  }
  
  # Apply removeBatchEffect to correct for known covariates
  known_covariates_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = covariates,
    design = design
  )
  
  # Add to the results list
  corrected_matrices$known_covariates <- known_covariates_matrix
  
  # ============================================================================
  # 3. Principal Component Regression (PCR)
  # ============================================================================
  message("Applying principal component regression...")
  
  # Filter genes with low variance
  gene_vars <- rowVars(normalized_matrix)
  high_var_genes <- which(gene_vars > quantile(gene_vars, 0.5))
  
  # Perform PCA on high-variance genes
  pca_data <- prcomp(t(normalized_matrix[high_var_genes, ]), scale = TRUE, center = TRUE)
  
  # Determine the number of significant PCs to remove
  # Create design matrix for the biological factor (group) only
  mod <- model.matrix(~ group, data = bulk_metadata)
  mod0 <- model.matrix(~ 1, data = bulk_metadata)
  
  # Estimate the number of surrogate variables
  n_sv <- num.sv(normalized_matrix, mod, method = "be", B = 20)
  message(paste("Estimated number of significant PCs:", n_sv))
  
  # Create a matrix of the PCs to regress out
  pcs_to_regress <- pca_data$x[, 1:n_sv, drop = FALSE]
  
  # Apply removeBatchEffect to remove the effect of significant PCs
  pc_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = pcs_to_regress,
    design = mod
  )
  
  # Add to the results list
  corrected_matrices$pc_correction <- pc_corrected_matrix
  
  # ============================================================================
  # 4. Surrogate Variable Analysis (SVA)
  # ============================================================================
  message("Applying surrogate variable analysis...")
  
  # Perform SVA
  sva_result <- sva(normalized_matrix, mod, mod0, n.sv = n_sv)
  
  # Extract surrogate variables
  sv <- sva_result$sv
  
  # Apply removeBatchEffect to remove the effect of surrogate variables
  sva_corrected_matrix <- limma::removeBatchEffect(
    normalized_matrix,
    covariates = sv,
    design = mod
  )
  
  # Add to the results list
  corrected_matrices$sva_correction <- sva_corrected_matrix
  
  # ============================================================================
  # 5. PEER (Probabilistic Estimation of Expression Residuals)
  # ============================================================================
  message("Applying PEER correction...")
  
  # Determine the number of PEER factors based on sample size
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
  message(paste("Using", n_factors, "PEER factors for", n_samples, "samples"))
  
  # Set up PEER model
  peer_model <- PEER()
  PEER_setPhenoMean(peer_model, t(normalized_matrix))
  
  # Add known covariates (group) to protect biological signal
  PEER_setNk(peer_model, n_factors)
  PEER_update(peer_model)
  
  # Get residuals (corrected data)
  residuals <- t(PEER_getResiduals(peer_model))
  rownames(residuals) <- rownames(normalized_matrix)
  colnames(residuals) <- colnames(normalized_matrix)
  
  # Add to the results list
  corrected_matrices$peer_correction <- residuals
  
  # ============================================================================
  # 6. RUVCorr
  # ============================================================================
  message("Applying RUVCorr correction...")
  
  # Find negative control genes (genes not in any module)
  module_cols <- grep("^in_module_", colnames(gene_metadata), value = TRUE)
  
  if (length(module_cols) == 0) {
    warning("No in_module columns found in gene_metadata; using all genes for RUVCorr")
    control_genes <- rep(TRUE, nrow(gene_metadata))
  } else {
    # Select genes that are not in any module
    control_genes <- rowSums(gene_metadata[, module_cols, drop = FALSE]) == 0
  }
  
  # If no control genes are found, use genes with lowest variance
  if (sum(control_genes) < 10) {
    warning("Fewer than 10 control genes found; using 10% of genes with lowest variance")
    control_genes <- gene_vars <= quantile(gene_vars, 0.1)
  }
  
  message(paste("Using", sum(control_genes), "negative control genes for RUVCorr"))
  
  # Set number of unwanted variation factors to remove (k)
  k <- min(5, ncol(bulk_counts) / 10)  # Use 5 or 10% of sample size, whichever is smaller
  
  # Apply RUVCorr
  ruv_result <- RUVcorr(Y = normalized_matrix, 
                        X = model.matrix(~ group, data = bulk_metadata),
                        ctl = control_genes, 
                        k = k)
  
  # Extract corrected matrix
  ruv_corrected_matrix <- ruv_result$Y_corr
  
  # Add to the results list
  corrected_matrices$ruv_corr <- ruv_corrected_matrix
  
  # ============================================================================
  # Return all corrected matrices
  # ============================================================================
  message("All correction methods applied successfully.")
  return(corrected_matrices)
}

# ============================================================================
# Example usage
# ============================================================================

#' Example of how to use the apply_all_corrections function
example_usage <- function() {
  # Load simulated data
  # Replace with your actual file path
  sim_data <- readRDS("path/to/your/simulated_data.rds")
  
  # Apply all correction methods
  corrected_data <- apply_all_corrections(sim_data)
  
  # Print dimensions of each corrected matrix
  for (method_name in names(corrected_data)) {
    matrix_dim <- dim(corrected_data[[method_name]])
    cat(sprintf("%s: %d genes × %d samples\n", 
                method_name, matrix_dim[1], matrix_dim[2]))
  }
  
  # Return the corrected data
  return(corrected_data)
}

# Uncomment to run the example
# results <- example_usage()