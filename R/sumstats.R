#' Simple wrapper for readr::read_tsv to use a gzipped file.
#'
#' @param filename A string ending in .tsv.gz or .tsv.bgz to be read.
read_gz_tsv <- function(filename){
  print(filename)
  readr::read_tsv(gzfile(filename))
}

#' Wrapper for read_gz_tsv to read all files matching pattern.
#'
#' @param file_pattern A string representing a glob pattern of files to match.
#' @param min_maf Float, minimum allele frequency to keep.
read_file_pattern <- function(file_pattern, min_maf = 0.001){
  files = Sys.glob(file_pattern)
  names = sapply(strsplit(sapply(
    strsplit(files, '/'), function(x){x[length(x)]}), "\\."), function(x){x[1]})
  files = stats::setNames(as.list(files), names)
  lapply(files, function(file){parse_ukbbss_neale(read_gz_tsv(file), min_maf)})
}

#' Reads single sumstats file in the format of the Neale lab UKBB analysis.
#'
#' Note that this converts beta values to a normalized scale if they are not
#' already on one. Note also that this reads the entire dataset into memory.
#' If the number of SNPs is
#' large, consider doing some preprocessing of your data to remove SNPs that
#' are unlikely to be used in the analysis.
#' Note also that this only keeps the normalized scale beta and SE due to memory
#' concerns. The number of samples and p-value can be backed out from these.
#'
#' @param sumstats Tibble. Must include ccolumns "minor_AF",
#'   "low_confidence_variant", "tstat", "n_completet_samples", "beta", "se"
#'   and "pval".
#' @param min_maf Float. Filter alleles with minor allele frequency less than
#'   min_maf.
parse_ukbbss_neale <- function(sumstats, min_maf = 0.001){
  # TODO(brielin): Use chunked version to filter on read.
  print("Parsing...")
  sumstats <- sumstats %>%
    dplyr::filter(minor_AF > min_maf, low_confidence_variant == FALSE) %>%
    dplyr::mutate(se_hat = 1/sqrt(n_complete_samples), beta_hat = tstat*se_hat) %>%
    dplyr::select(variant, beta_hat, se_hat)
  return(sumstats)
}

#' Parses a list of sumstats tibbles, returning a list of matrices.
#'
#' This parses a list of sumstats tibbles into a list of matrices containing
#' the relevant data for future analysis, ie that required by `fit_sumstats`.
#'
#' @param sumstats A list where each entry is a tibble of sumstats.
parse_sumstats_multi_trait <- function(sumstats){
  sumstats <- dplyr::bind_rows(sumstats, .id = "trait")
  beta_hat <- sumstats %>% dplyr::select("variant", "trait", "beta_hat") %>%
    tidyr::pivot_wider(names_from = "trait", values_from = "beta_hat", id_cols = "variant") %>%
    tibble::column_to_rownames("variant")
  se_hat <- sumstats %>% dplyr::select("variant", "trait", "se_hat") %>%
    tidyr::pivot_wider(names_from = "trait", values_from = "se_hat", id_cols = "variant") %>%
    tibble::column_to_rownames("variant")
  list("beta_hat" = beta_hat, "se_hat" = se_hat)
}

#' Reads a pattern of summary statistics formatted according to the Neale lab.
#'
#' @param file_pattern A glob string specifying the files to read.
#' @param min_maf Float, minimum allele frequency to keep.
read_ukbbss_neale <- function(file_pattern, min_maf = 0.001){
  sumstats <- read_file_pattern(file_pattern, min_maf)
  parse_sumstats_multi_trait(sumstats)
}
