# ==============================================================================
# ATLANTIC HERRING
# Filter genomic data and estimate pairwise FST
# ==============================================================================
#
# Input:
#   60.Neff.allele_freq
#
# Pipeline:
#   1. Read allele-frequency data.
#   2. Exclude non-target pools listed below.
#   3. Retain chromosomes chr1-chr26.
#   4. Select approximately neutral SNPs using:
#        - MAF proxy > 0.01
#        - 0.02 <= SD(allele_freq) <= 0.08
#        - data available for >= 90% of populations
#   5. Thin to one SNP per 1-kb window.
#   6. Estimate pairwise FST among the 53 retained samples.
#   7. Save pairwise matrices, long-format results, and filtering summaries.
#
# Multilocus FST estimator:
#
#                      sum(Ht - Hs)
#                FST = ------------
#                         sum(Ht)
#
# where:
#   H1   = 2 * p1 * (1 - p1)
#   H2   = 2 * p2 * (1 - p2)
#   Hs   = (H1 + H2) / 2
#   pbar = (p1 + p2) / 2
#   Ht   = 2 * pbar * (1 - pbar)
#
# ==============================================================================


# ==============================================================================
# 0. SETUP AND PACKAGES
# ==============================================================================

library(data.table)
library(matrixStats)

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

input_file <- "60.Neff.allele_freq"

# Minimum proportion of populations with data for each SNP.
min_prop_pop <- 0.90

# Criteria for SNPs with low among-population differentiation.
sd_min <- 0.02
sd_max <- 0.08

# Minimum minor-allele-frequency proxy.
maf_min <- 0.01

# Main chromosomes.
chr_keep <- paste0("chr", 1:26)

# Thinning window size in base pairs.
window_size_bp <- 1000L


# ==============================================================================
# 2. READ ALLELE-FREQUENCY DATA
# ==============================================================================

allele_freq <- fread(input_file)


# ==============================================================================
# 3. IDENTIFY POPULATION COLUMNS
# ==============================================================================

population_cols <- setdiff(
  names(allele_freq),
  c("CHROM", "POS")
)


# ==============================================================================
# 4. EXCLUDE NON-TARGET POOLS
# ==============================================================================

non_target_cols <- c(
  "HWS1_Japan_SeaOfJapan",
  "HWS2_PechoraSea_BarentsSea",
  "HWS3_WhiteSea_WhiteSea",
  "HWS4_KandalakshaBay_WhiteSea",
  "HWS5_KandalakshaBay_WhiteSea",
  "HWS6_Balsfjord_Atlantic",
  "PB8_Pacific_Pacific_Spring"
)

stopifnot(all(non_target_cols %in% population_cols))

atlantic_population_cols <- setdiff(population_cols, non_target_cols)

stopifnot(length(atlantic_population_cols) == 53)


# ==============================================================================
# 5. MINIMUM NUMBER OF POPULATIONS WITH DATA
# ==============================================================================

min_n_pop <- ceiling(length(atlantic_population_cols) * min_prop_pop)


# ==============================================================================
# 6. FILTER SNPs AND APPLY 1-kb THINNING
# ==============================================================================
#
# Filtering is performed chromosome by chromosome to avoid creating a
# ~15-million x 53 matrix in memory.
#
# ==============================================================================

neutral_snps_list <- vector("list", length(chr_keep))
filter_summary_list <- vector("list", length(chr_keep))

for (k in seq_along(chr_keep)) {
  chr <- chr_keep[k]
  
  # ---------------------------------------------------------------------------
  # Extract chromosome
  # ---------------------------------------------------------------------------
  
  chromosome_data <- allele_freq[
    CHROM == chr,
    c("CHROM", "POS", atlantic_population_cols),
    with = FALSE
  ]
  
  n_original <- nrow(chromosome_data)
  
  # ---------------------------------------------------------------------------
  # Convert only the current chromosome to a matrix
  # ---------------------------------------------------------------------------
  
  chromosome_freq_matrix <- as.matrix(
    chromosome_data[, ..atlantic_population_cols]
  )
  
  # Number of populations with data.
  n_pop <- matrixStats::rowCounts(!is.na(chromosome_freq_matrix))
  
  # Mean allele frequency across populations.
  mean_p <- matrixStats::rowMeans2(chromosome_freq_matrix, na.rm = TRUE)
  
  # MAF proxy based on pool-level allele frequencies:
  #   maf_proxy = min(mean_p, 1 - mean_p)
  maf_proxy <- pmin(mean_p, 1 - mean_p)
  
  # Standard deviation of allele frequency across populations.
  sd_pop <- matrixStats::rowSds(chromosome_freq_matrix, na.rm = TRUE)
  
  # ---------------------------------------------------------------------------
  # Select approximately neutral loci
  # ---------------------------------------------------------------------------
  
  keep <- (
    n_pop >= min_n_pop &
      is.finite(maf_proxy) &
      is.finite(sd_pop) &
      maf_proxy > maf_min &
      sd_pop >= sd_min &
      sd_pop <= sd_max
  )
  
  n_filtered <- sum(keep)
  
  retained_snps <- chromosome_data[keep]
  
  retained_snps[, n_pop := n_pop[keep]]
  retained_snps[, maf_proxy := maf_proxy[keep]]
  retained_snps[, sd_pop := sd_pop[keep]]
  
  # Release the large chromosome-level matrix before thinning.
  rm(
    chromosome_freq_matrix,
    n_pop,
    mean_p,
    maf_proxy,
    sd_pop,
    keep,
    chromosome_data
  )
  invisible(gc())
  # ---------------------------------------------------------------------------
  # Thin to one eligible SNP per 1-kb window
  # ---------------------------------------------------------------------------
  
  if (nrow(retained_snps) > 0) {
    setorder(retained_snps, POS)
    
    # Windows are 1-1000, 1001-2000, 2001-3000, ...
    retained_snps[, bin_1kb := (POS - 1L) %/% window_size_bp]
    
    # Retain the first eligible SNP in each window.
    thinned_snps <- retained_snps[!duplicated(bin_1kb)]
    thinned_snps[, bin_1kb := NULL]
  } else {
    thinned_snps <- retained_snps
  }
  
  n_thinned <- nrow(thinned_snps)
  
  # Store chromosome-level results.
  neutral_snps_list[[k]] <- thinned_snps
  
  filter_summary_list[[k]] <- data.table(
    CHROM = chr,
    n_original = n_original,
    n_neutral = n_filtered,
    n_after_1kb_thinning = n_thinned
  )
  
  rm(retained_snps, thinned_snps)
  invisible(gc())
}


# ==============================================================================
# 7. COMBINE RESULTS ACROSS CHROMOSOMES
# ==============================================================================

neutral_snps <- rbindlist(
  neutral_snps_list,
  use.names = TRUE
)

filter_summary <- rbindlist(filter_summary_list)


# ==============================================================================
# 8. SUMMARIZE RETAINED SNPs
# ==============================================================================


# ==============================================================================
# 9. CREATE SNP x POPULATION MATRIX
# ==============================================================================
#
# At this stage the matrix is manageable because only filtered and thinned
# SNPs are retained.
#
# ==============================================================================

neutral_freq_matrix <- as.matrix(neutral_snps[, ..atlantic_population_cols])


# ==============================================================================
# 10. INITIALIZE PAIRWISE MATRICES
# ==============================================================================

n_populations <- ncol(neutral_freq_matrix)

fst_pairwise <- matrix(
  NA_real_,
  nrow = n_populations,
  ncol = n_populations,
  dimnames = list(atlantic_population_cols, atlantic_population_cols)
)

n_loci_pairwise <- matrix(
  0L,
  nrow = n_populations,
  ncol = n_populations,
  dimnames = list(atlantic_population_cols, atlantic_population_cols)
)


# ==============================================================================
# 11. ESTIMATE PAIRWISE FST
# ==============================================================================

for (i in seq_len(n_populations - 1L)) {
  for (j in (i + 1L):n_populations) {
    # Allele frequencies in the two populations.
    p1 <- neutral_freq_matrix[, i]
    p2 <- neutral_freq_matrix[, j]
    
    # Use only loci with finite allele frequencies in both populations.
    ok <- is.finite(p1) & is.finite(p2)
    
    p1 <- p1[ok]
    p2 <- p2[ok]
    
    n_loci <- length(p1)
    
    n_loci_pairwise[i, j] <- n_loci
    n_loci_pairwise[j, i] <- n_loci
    
    if (n_loci == 0L) {
      next
    }
    
    # Expected heterozygosity within each population.
    H1 <- 2 * p1 * (1 - p1)
    H2 <- 2 * p2 * (1 - p2)
    
    # Average heterozygosity within populations.
    Hs <- (H1 + H2) / 2
    
    # Mean allele frequency across the pair.
    pbar <- (p1 + p2) / 2
    
    # Expected heterozygosity in the pooled pair.
    Ht <- 2 * pbar * (1 - pbar)
    
    # Multilocus FST:
    #
    #                    sum(Ht - Hs)
    #              FST = ------------
    #                       sum(Ht)
    #
    # This is a ratio of sums, not the arithmetic mean of locus-specific FST.
    denominator <- sum(Ht, na.rm = TRUE)
    
    if (is.finite(denominator) && denominator > 0) {
      fst <- sum(Ht - Hs, na.rm = TRUE) / denominator
    } else {
      fst <- NA_real_
    }
    
    # Store symmetrically.
    fst_pairwise[i, j] <- fst
    fst_pairwise[j, i] <- fst
  }
}


# ==============================================================================
# 12. SET FST MATRIX DIAGONAL
# ==============================================================================

diag(fst_pairwise) <- 0


# ==============================================================================
# 13. CONVERT PAIRWISE FST TO LONG FORMAT
# ==============================================================================
#
# Keep one observation per pair because the FST matrix is symmetric.
#
# ==============================================================================

pair_indices <- which(
  upper.tri(fst_pairwise),
  arr.ind = TRUE
)

fst_long <- data.table(
  pop1 = rownames(fst_pairwise)[pair_indices[, 1]],
  pop2 = colnames(fst_pairwise)[pair_indices[, 2]],
  FST = fst_pairwise[pair_indices],
  n_snps = n_loci_pairwise[pair_indices]
)

# Sort from highest to lowest differentiation.
setorder(fst_long, -FST)


# ==============================================================================
# 14. SAVE OUTPUTS
# ==============================================================================

fst_out <- data.table(
  population = rownames(fst_pairwise)
)

fst_out <- cbind(
  fst_out,
  as.data.table(fst_pairwise)
)

fwrite(
  fst_out,
  "Herring_FST_pairwise_neutral.csv"
)

fwrite(
  fst_long,
  "Herring_FST_pairwise_neutral_long.csv"
)

fwrite(
  filter_summary,
  "Herring_neutral_SNP_filter_summary.csv"
)

saveRDS(
  list(
    fst_pairwise = fst_pairwise,
    fst_long = fst_long,
    n_loci_pairwise = n_loci_pairwise,
    neutral_snps = neutral_snps,
    filter_summary = filter_summary,
    atlantic_population_cols = atlantic_population_cols,
    parameters = list(
      min_prop_pop = min_prop_pop,
      sd_min = sd_min,
      sd_max = sd_max,
      maf_min = maf_min
    )
  ),
  file = "Herring_FST_neutral_analysis.rds"
)