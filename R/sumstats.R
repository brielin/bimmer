#' Simple wrapper for readr::read_tsv to use a gzipped file.
#'
#' @param filename A string ending in .tsv.gz or .tsv.bgz to be read.
read_gz_tsv <- function(filename) {
  print(filename)
  return(readr::read_tsv(gzfile(filename)))
}

#' Wrapper for read_gz_tsv to read all files matching pattern.
#'
#' One of file_pattern and file_list must be provided.
#'
#' @param file_pattern A string representing a glob pattern of files to match.
#' @param file_list A list of files to open.
read_files <- function(file_pattern = NULL, file_list = NULL) {
  if(!is.null(file_pattern)){
    files <- Sys.glob(file_pattern)
  } else if (!is.null(file_list)){
    files <- file_list
  }
  else {
    stop("One of file_pattern or file_list must be provided!")
  }

  names <- sapply(strsplit(sapply(
    strsplit(files, "/"), function(x) {
      x[length(x)]
    }
  ), "\\."), function(x) {
    x[1]
  })
  files <- stats::setNames(as.list(files), names)
  return(lapply(files, function(file) {
    parse_ukbbss_neale(read_gz_tsv(file))
  }))
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
parse_ukbbss_neale <- function(sumstats) {
  # TODO(brielin): Use chunked version to filter on read.
  print("Parsing...")
  # Manually set low confidence to NA instead of dropping them to keep the SNP
  # lists identical.
  sumstats$tstat[sumstats$low_confidence_variant == TRUE] <- NA
  sumstats <- sumstats %>%
    dplyr::filter(!duplicated(sumstats$rsid)) %>%
    dplyr::mutate(
      se_hat = 1 / sqrt(n_complete_samples), beta_hat = tstat * se_hat) %>%
    dplyr::select(rsid, beta_hat, se_hat)
  return(sumstats)
}

#' Parses a list of sumstats tibbles, returning a list of matrices.
#'
#' This parses a list of sumstats tibbles into a list of matrices containing
#' the relevant data for future analysis, ie that required by `fit_sumstats`.
#'
#' @param sumstats A list where each entry is a tibble of sumstats.
parse_sumstats_multi_trait <- function(sumstats) {
  sumstats <- dplyr::bind_rows(sumstats, .id = "trait")
  beta_hat <- sumstats %>%
    dplyr::select("rsid", "trait", "beta_hat") %>%
    tidyr::pivot_wider(
      names_from = "trait", values_from = "beta_hat", id_cols = "rsid") %>%
    tibble::column_to_rownames("rsid")
  se_hat <- sumstats %>%
    dplyr::select("rsid", "trait", "se_hat") %>%
    tidyr::pivot_wider(
      names_from = "trait", values_from = "se_hat", id_cols = "rsid") %>%
    tibble::column_to_rownames("rsid")
  return(list("beta_hat" = beta_hat, "se_hat" = se_hat))
}

#' Reads a pattern of summary statistics formatted according to the Neale lab.
#'
#' @param file_pattern A glob string specifying the files to read.
#' @param chunk_size Number of phenotypes to pivot at a time.
read_ukbbss_neale <- function(file_pattern = NULL, file_list = NULL,
                              chunk_size = 20) {
  sumstats <- read_files(file_pattern, file_list)

  # Ideally, we would rbind sumstats and then pivot off the columns we want
  # to get matrices. However they are too big and this breaks tidyr. Instead
  # we want to pop off chunks of sumstats at a time and pivot them, without
  # having to store multiple copies of the data. Since R protects objects in
  # functions via the namespace, we need to modify the sumstats object in
  # the place that it was created.
  pivot_chunk <- function(nread) {
    print("Pivot chunk!")
    ss_chunk <- sumstats[1:nread]
    sumstats[1:nread] <<- NULL # Remove the elements we just popped.
    print(length(sumstats))
    ss_chunk <- dplyr::bind_rows(ss_chunk, .id = "trait")
    beta_chunk <- ss_chunk %>%
      dplyr::select("rsid", "trait", "beta_hat") %>%
      tidyr::pivot_wider(
        names_from = "trait", values_from = "beta_hat", id_cols = "rsid")
    se_chunk <- ss_chunk %>%
      dplyr::select("rsid", "trait", "se_hat") %>%
      tidyr::pivot_wider(
        names_from = "trait", values_from = "se_hat", id_cols = "rsid")
    list("beta" = beta_chunk, "se" = se_chunk)
  }
  num_stats <- length(sumstats)
  if (num_stats %% chunk_size > 0) {
    nread_seq <- c(rep(chunk_size, num_stats %/% chunk_size),
                   num_stats %% chunk_size)
  } else {
    nread_seq <- rep(chunk_size, num_stats %/% chunk_size)
  }
  print(length(sumstats))
  pivoted_chunks <- purrr::map(nread_seq, pivot_chunk)
  pivoted_chunks <- purrr::transpose(pivoted_chunks)
  # This is a hack to check that the rsids in each chunk match.
  if (all(sapply(pivoted_chunks$beta, function(x) {
    all.equal(x$rsid, pivoted_chunks$beta[[1]]$rsid)
  }))) {
    return(list(
      "beta_hat" = dplyr::bind_cols(pivoted_chunks$beta) %>%
        tibble::column_to_rownames("rsid") %>%
        dplyr::select(-dplyr::starts_with("rsid")),
      "se_hat" = dplyr::bind_cols(pivoted_chunks$se) %>%
        tibble::column_to_rownames("rsid") %>%
        dplyr::select(-dplyr::starts_with("rsid"))
    ))
  } else {
    warning("The variant IDs in each chunk do not match! Performing a full join.
            This may take a while.")
    return(list(
      "beta_hat" = plyr::join_all(
        pivoted_chunks$beta, by = "rsid", type = "full", match = "first") %>%
        tibble::column_to_rownames("rsid"),
      "se_hat" = plyr::join_all(
        pivoted_chunks$se, by = "rsid", type = "full", match = "first") %>%
        tibble::column_to_rownames("rsid")
    ))
  }
}

#' Reads a file pattern corresponding to files with potential instruments.
#'
#' @param file_pattern A file pattern corresponding to a set of files, one for
#'   each phenotype. Each file should contain a list of potential SNPs to use as
#'   instruments for that phenotype.
read_snp_list <- function(file_pattern) {
  files <- Sys.glob(file_pattern)
  names <- sapply(strsplit(sapply(
    strsplit(files, "/"), function(x) {
      x[length(x) - 1]
    }
  ), "\\."), function(x) {
    x[1]
  })
  files <- stats::setNames(as.list(files), names)
  return(lapply(files, function(file) {
    readLines(file)
  }))
}
