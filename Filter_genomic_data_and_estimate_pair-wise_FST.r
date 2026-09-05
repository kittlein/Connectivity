# Filter genomic data and estimate pair-wise FST for Atlantic herring.
# Retains chr1-chr26, filters loci, thins SNPs, and estimates multilocus FST.

library(data.table)
library(matrixStats)

# Define filtering criteria.
# SNPs frecuencies from https://doi.org/10.5061/dryad.pnvx0k6kr
input_file <- "60.Neff.freq"

min_prop_pop <- 0.90
sd_min <- 0.02
sd_max <- 0.08
maf_min <- 0.01
chr_keep <- paste0("chr", 1:26)
window_size_bp <- 1000L

# Read allele-frequency data.
allele_freq <- fread(input_file)

# Identify population columns.
population_cols <- setdiff(
  names(allele_freq),
  c("CHROM", "POS")
)

# Exclude non-target pools.
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

atlantic_population_cols <- setdiff(
  population_cols,
  non_target_cols
)

stopifnot(length(atlantic_population_cols) == 53)

# Require data for at least 90% of retained populations.
min_n_pop <- ceiling(
  length(atlantic_population_cols) * min_prop_pop
)

# Filter each chromosome separately to limit memory use.
neutral_snps_list <- vector("list", length(chr_keep))
filter_summary_list <- vector("list", length(chr_keep))

for (k in seq_along(chr_keep)) {
  chr <- chr_keep[k]

  chromosome_data <- allele_freq[
    CHROM == chr,
    c("CHROM", "POS", atlantic_population_cols),
    with = FALSE
  ]

  n_original <- nrow(chromosome_data)

  chromosome_freq_matrix <- as.matrix(
    chromosome_data[, ..atlantic_population_cols]
  )

  # Calculate completeness, mean frequency, MAF proxy, and frequency SD.
  n_pop <- matrixStats::rowCounts(
    !is.na(chromosome_freq_matrix)
  )

  mean_p <- matrixStats::rowMeans2(
    chromosome_freq_matrix,
    na.rm = TRUE
  )

  maf_proxy <- pmin(
    mean_p,
    1 - mean_p
  )

  sd_pop <- matrixStats::rowSds(
    chromosome_freq_matrix,
    na.rm = TRUE
  )

  # Retain loci meeting completeness, MAF, and differentiation criteria.
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

  # Thin to the first eligible SNP in each 1-kb genomic window.
  if (nrow(retained_snps) > 0) {
    setorder(retained_snps, POS)

    retained_snps[
      ,
      bin_1kb := (POS - 1L) %/% window_size_bp
    ]

    thinned_snps <- retained_snps[
      !duplicated(bin_1kb)
    ]

    thinned_snps[
      ,
      bin_1kb := NULL
    ]
  } else {
    thinned_snps <- retained_snps
  }

  n_thinned <- nrow(thinned_snps)

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

# Combine retained SNPs and chromosome-level filtering summaries.
neutral_snps <- rbindlist(
  neutral_snps_list,
  use.names = TRUE
)

filter_summary <- rbindlist(
  filter_summary_list
)

# Create the final SNP-by-population frequency matrix.
neutral_freq_matrix <- as.matrix(
  neutral_snps[, ..atlantic_population_cols]
)

# Initialize pair-wise FST and loci-count matrices.
n_populations <- ncol(neutral_freq_matrix)

fst_pairwise <- matrix(
  NA_real_,
  nrow = n_populations,
  ncol = n_populations,
  dimnames = list(
    atlantic_population_cols,
    atlantic_population_cols
  )
)

n_loci_pairwise <- matrix(
  0L,
  nrow = n_populations,
  ncol = n_populations,
  dimnames = list(
    atlantic_population_cols,
    atlantic_population_cols
  )
)

# Estimate multilocus FST as sum(Ht - Hs) / sum(Ht).
for (i in seq_len(n_populations - 1L)) {
  for (j in (i + 1L):n_populations) {
    p1 <- neutral_freq_matrix[, i]
    p2 <- neutral_freq_matrix[, j]

    ok <- is.finite(p1) & is.finite(p2)

    p1 <- p1[ok]
    p2 <- p2[ok]

    n_loci <- length(p1)

    n_loci_pairwise[i, j] <- n_loci
    n_loci_pairwise[j, i] <- n_loci

    if (n_loci == 0L) {
      next
    }

    H1 <- 2 * p1 * (1 - p1)
    H2 <- 2 * p2 * (1 - p2)
    Hs <- (H1 + H2) / 2

    pbar <- (p1 + p2) / 2
    Ht <- 2 * pbar * (1 - pbar)

    denominator <- sum(
      Ht,
      na.rm = TRUE
    )

    if (is.finite(denominator) && denominator > 0) {
      fst <- sum(
        Ht - Hs,
        na.rm = TRUE
      ) / denominator
    } else {
      fst <- NA_real_
    }

    fst_pairwise[i, j] <- fst
    fst_pairwise[j, i] <- fst
  }
}

diag(fst_pairwise) <- 0

# extract samples east of 20W
localities = fread("SNP_samples_herring.csv")
iloc <- localities$sample[localities$sample %in% colnames(fst_pairwise)]

# write FST data.table
fwrite(data.table(fst_pairwise[iloc, iloc]), "paired_Fst.csv")
