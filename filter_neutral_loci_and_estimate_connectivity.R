## ------------------------------------------------------------
## Export herring SNP genotypes from FinePop2 to GENEPOP
## retaining only putatively neutral loci
## ------------------------------------------------------------

library(FinePop2)

data("herring")

## herring$genepop is already GENEPOP-format text
gp <- herring$genepop

## Make sure it is a vector of lines
if (length(gp) == 1) {
  gp <- unlist(strsplit(gp, "\n", fixed = TRUE))
}

gp <- gsub("\r", "", gp)
gp <- gp[nchar(trimws(gp)) > 0]

## ------------------------------------------------------------
## Function to subset loci in a GENEPOP text object
## ------------------------------------------------------------

subset_genepop_loci <- function(gp_lines, loci_to_remove, outfile) {
  
  pop_lines <- grep("^pop$", trimws(gp_lines), ignore.case = TRUE)
  if (length(pop_lines) == 0) stop("No POP lines found in the GENEPOP file.")
  
  first_pop <- pop_lines[1]
  
  title <- gp_lines[1]
  locus_block <- gp_lines[2:(first_pop - 1)]
  
  ## FinePop2 herring usually has one locus per line,
  ## but this also handles comma-separated locus headers.
  if (any(grepl(",", locus_block))) {
    loci <- trimws(unlist(strsplit(paste(locus_block, collapse = ","), ",")))
  } else {
    loci <- trimws(locus_block)
  }
  
  ## Normalize names to match both "1025.1-149" and "Cha_1025.1-149"
  normalize_locus <- function(x) {
    x <- trimws(x)
    x <- gsub("^Cha_", "", x)
    x <- gsub("\\s*\\(.*\\)$", "", x)
    x
  }
  
  loci_norm <- normalize_locus(loci)
  remove_norm <- normalize_locus(loci_to_remove)
  
  keep <- !(loci_norm %in% remove_norm)
  
  message("Original number of loci: ", length(loci))
  message("Loci removed: ", sum(!keep))
  message("Loci retained: ", sum(keep))
  
  if (sum(!keep) != length(unique(remove_norm))) {
    warning(
      "Some loci listed for removal were not found in the GENEPOP header. ",
      "Check locus names with setdiff(remove_norm, loci_norm)."
    )
    print(setdiff(remove_norm, loci_norm))
  }
  
  data_block <- gp_lines[first_pop:length(gp_lines)]
  
  new_data_block <- character(length(data_block))
  
  for (i in seq_along(data_block)) {
    
    line <- data_block[i]
    
    if (grepl("^pop$", trimws(line), ignore.case = TRUE)) {
      new_data_block[i] <- "POP"
      next
    }
    
    ## Individual lines: individual_name, genotype1 genotype2 ...
    parts <- strsplit(line, ",", fixed = TRUE)[[1]]
    
    if (length(parts) < 2) {
      stop("Unexpected individual line without comma: ", line)
    }
    
    indiv_id <- trimws(parts[1])
    geno_string <- trimws(paste(parts[-1], collapse = ","))
    genotypes <- unlist(strsplit(geno_string, "\\s+"))
    
    if (length(genotypes) != length(loci)) {
      stop(
        "Number of genotypes does not match number of loci for individual ",
        indiv_id, ". Found ", length(genotypes),
        " genotypes and ", length(loci), " loci."
      )
    }
    
    new_data_block[i] <- paste0(
      indiv_id, ", ",
      paste(genotypes[keep], collapse = " ")
    )
  }
  
  new_gp <- c(
    title,
    loci[keep],
    new_data_block
  )
  
  writeLines(new_gp, con = outfile)
  
  invisible(
    list(
      outfile = outfile,
      original_loci = loci,
      retained_loci = loci[keep],
      removed_loci = loci[!keep]
    )
  )
}

global_outliers_limborg2012 <-  c(
  "Cha_10193.1-449",
  "Cha_1068.2-349",
  "Cha_1165.2-123",
  "Cha_13178.2-124",
  "Cha_13259.1-167",
  "Cha_143.1-185",
  "Cha_15105.2-341",
  "Cha_15360.2-279",
  "Cha_15984.1-275",
  "Cha_16330.7-357",
  "Cha_2814.1-396",
  "Cha_2884.1-367",
  "Cha_297.1-93",
  "Cha_381.2-437",
  "Cha_5534.1-506",
  "Cha_5541.1-273",
  "Cha_11896.1-201",
  "Cha_13376.1-166",
  "Cha_16060.1-279",
  "Cha_462.3-102",
  "Cha_10361.2-383",
  "Cha_1143.1-484",
  "Cha_12771.1-298",
  "Cha_13097.1-122",
  "Cha_15034.1-201",
  "Cha_15659.7-503",
  "Cha_15898.2-568",
  "Cha_160.1-805",
  "Cha_535.2-394",
  "Cha_5625.1-135",
  "Cha_9634.1-256",
  "Cha_1400.3-301",
  "Cha_3888.1-826",
  "Cha_688.1-238",
  "Cha_693.2-263",
  "Cha_15964.1-332",
  "Cha_318.1-301",
  "Cha_7456.1-168",
  "Cha_8760.1-243"
)

outliers_strict <- c(global_outliers_limborg2012,
                     "Cha_10193.1-449",
                     "Cha_1068.2-349",
                     "Cha_1165.2-123",
                     "Cha_13178.2-124",
                     "Cha_13259.1-167",
                     "Cha_143.1-185",
                     "Cha_15105.2-341",
                     "Cha_15360.2-279",
                     "Cha_15984.1-275",
                     "Cha_16330.7-357",
                     "Cha_2814.1-396",
                     "Cha_2884.1-367",
                     "Cha_297.1-93",
                     "Cha_381.2-437",
                     "Cha_5534.1-506",
                     "Cha_5541.1-273",
                     "Cha_11896.1-201",
                     "Cha_13376.1-166",
                     "Cha_16060.1-279",
                     "Cha_462.3-102",
                     "Cha_10361.2-383",
                     "Cha_1143.1-484",
                     "Cha_12771.1-298",
                     "Cha_13097.1-122",
                     "Cha_15034.1-201",
                     "Cha_15659.7-503",
                     "Cha_15898.2-568",
                     "Cha_160.1-805",
                     "Cha_535.2-394",
                     "Cha_5625.1-135",
                     "Cha_9634.1-256",
                     "Cha_1400.3-301",
                     "Cha_3888.1-826",
                     "Cha_688.1-238",
                     "Cha_693.2-263",
                     "Cha_15964.1-332",
                     "Cha_318.1-301",
                     "Cha_7456.1-168",
                     "Cha_8760.1-243"
)

res_neutral_242 <- subset_genepop_loci(
  gp_lines = gp,
  loci_to_remove = outliers_strict,
  outfile = "herring_neutral_242_Limborg2012.genepop"
)

library(diveRsity)

# GENEPOP file containing the neutral loci
infile <- "herring_neutral_242_Limborg2012.genepop"

mig <- divMigrate(
  infile = infile,
  stat = "all",           # computes several metrics if supported by your version
  boots = 1000,           # bootstrap replicates for confidence intervals
  plot_network = FALSE,   # automatically plots the network
  filter_threshold = 0,   # displays all connections
  para = TRUE
)

# Copy of the matrix
con2 <- as.matrix(mig$dRelMig)

saveRDS(con2, "conectNeutralLoci.rds")
