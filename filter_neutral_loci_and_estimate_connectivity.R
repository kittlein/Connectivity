# Curated workflow for neutral-locus filtering and directional connectivity estimation.

library(FinePop2)
library(diveRsity)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

output_genepop_file <- "herring_neutral_265.gen"
output_locus_file <- "neutral_locus_names.csv"
output_connectivity_matrix_file <- "neutral_connectivity_matrix.rds"
output_connectivity_edges_file <- "neutral_connectivity_edges.csv"

# Global outlier loci reported in Limborg et al. (2012), Table 2.
global_outlier_loci <- c(
  "Cha_1025.1-149",
  "Cha_10733.1-102",
  "Cha_1170.1-250",
  "Cha_13197.4-115",
  "Cha_13371.3-81",
  "Cha_14331.1-140",
  "Cha_1513.1-91",
  "Cha_15389.3-101",
  "Cha_15984.1-275",
  "Cha_16330.7-357",
  "Cha_2814.1-396",
  "Cha_2884.1-367",
  "Cha_297.1-93",
  "Cha_381.2-437",
  "Cha_5534.1-506",
  "Cha_5541.1-273"
)

standardize_locus_name <- function(locus_name) {
  locus_name |>
    str_remove("\\s*\\(.*?\\)") |>
    str_trim() |>
    str_remove("^Cha_") |>
    str_remove("^Ch_")
}

load_herring_genepop <- function() {
  data("herring", package = "FinePop2", envir = environment())

  if (!exists("herring", inherits = FALSE)) {
    stop("The `herring` dataset was not found in the FinePop2 package.")
  }

  genepop_lines <- get("herring", inherits = FALSE)

  if (length(genepop_lines) == 1) {
    genepop_lines <- strsplit(genepop_lines, "\\r?\\n")[[1]]
  }

  genepop_lines[nzchar(genepop_lines)]
}

parse_genepop_structure <- function(genepop_lines) {
  pop_index <- which(str_detect(genepop_lines, regex("^pop$", ignore_case = TRUE)))[1]

  if (is.na(pop_index)) {
    stop("Could not identify population blocks in the GENEPOP text.")
  }

  locus_lines <- genepop_lines[2:(pop_index - 1)]
  locus_names <- locus_lines |>
    paste(collapse = " ") |>
    str_split(",") |>
    unlist() |>
    str_trim() |>
    discard(~ .x == "")

  list(
    header = genepop_lines[1],
    locus_names = locus_names,
    population_lines = genepop_lines[pop_index:length(genepop_lines)]
  )
}

split_genotype_string <- function(genotype_string, n_loci) {
  trimmed_string <- str_squish(genotype_string)
  genotype_tokens <- str_split(trimmed_string, "\\s+")[[1]]
  genotype_tokens <- genotype_tokens[genotype_tokens != ""]

  if (length(genotype_tokens) == n_loci) {
    return(genotype_tokens)
  }

  compact_string <- str_remove_all(trimmed_string, "\\s+")

  if (nchar(compact_string) %% n_loci != 0) {
    stop("Genotype string could not be split into locus-wise values.")
  }

  token_width <- nchar(compact_string) / n_loci
  substring(
    compact_string,
    first = seq(1, nchar(compact_string), by = token_width),
    last = seq(token_width, nchar(compact_string), by = token_width)
  )
}

filter_genepop_loci <- function(genepop_lines, loci_to_remove) {
  genepop_parsed <- parse_genepop_structure(genepop_lines)

  standardized_loci <- standardize_locus_name(genepop_parsed$locus_names)
  keep_index <- !standardized_loci %in% standardize_locus_name(loci_to_remove)

  retained_loci <- genepop_parsed$locus_names[keep_index]
  retained_lines <- map_chr(genepop_parsed$population_lines, function(line) {
    if (str_detect(line, regex("^pop$", ignore_case = TRUE))) {
      return("Pop")
    }

    line_parts <- str_split_fixed(line, ",", n = 2)
    sample_id <- str_trim(line_parts[1])
    genotype_string <- str_trim(line_parts[2])
    genotype_tokens <- split_genotype_string(genotype_string, length(genepop_parsed$locus_names))
    retained_tokens <- genotype_tokens[keep_index]

    paste0(sample_id, ", ", paste(retained_tokens, collapse = " "))
  })

  filtered_lines <- c(
    genepop_parsed$header,
    paste(retained_loci, collapse = ", "),
    retained_lines
  )

  list(
    genepop_lines = filtered_lines,
    retained_loci = retained_loci,
    removed_loci = genepop_parsed$locus_names[!keep_index]
  )
}

estimate_directional_connectivity <- function(genepop_file) {
  divMigrate(
    infile = genepop_file,
    stat = "gst",
    filter_threshold = 0,
    boots = 1000,
    plot_network = FALSE
  )
}

herring_genepop_lines <- load_herring_genepop()
filtered_genepop <- filter_genepop_loci(herring_genepop_lines, global_outlier_loci)

writeLines(filtered_genepop$genepop_lines, output_genepop_file)
write.csv(
  data.frame(locus_name = filtered_genepop$retained_loci),
  output_locus_file,
  row.names = FALSE
)

connectivity_results <- estimate_directional_connectivity(output_genepop_file)
neutral_connectivity_matrix <- connectivity_results$gRelMig

diag(neutral_connectivity_matrix) <- 0
saveRDS(neutral_connectivity_matrix, output_connectivity_matrix_file)

connectivity_edges <- as.data.frame(as.table(neutral_connectivity_matrix)) |>
  rename(from_population = Var1, to_population = Var2, relative_migration = Freq) |>
  filter(from_population != to_population)

write.csv(connectivity_edges, output_connectivity_edges_file, row.names = FALSE)

message("Saved filtered GENEPOP file to: ", output_genepop_file)
message("Saved retained loci to: ", output_locus_file)
message("Saved connectivity matrix to: ", output_connectivity_matrix_file)
message("Saved directional connectivity edges to: ", output_connectivity_edges_file)