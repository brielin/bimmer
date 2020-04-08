#!/bin/bash
#
# This script downloads the separated male and female UK biobank summary
# statistics from the Neale lab dropbox by downloading every specified
# UKBB field from a TSV where the first column in UKBB field number. It
# Also downloads all FinnGen phenotype codes.
#
# Args:
#   $1: Manifest location.
#   $2: TSV specifying the UKBB fields to use.
#   $3: Desired download directory.
#   $4: Maximum number of files to download.
#   $5: File prefix to write male/female file lists to.



MANIFEST=$1
SELECTED_FIELDS=$2
DOWNLOAD_DIR=$3
MAX_DOWNLOADS=$4

[ ! -d "${DOWNLOAD_DIR}" ] && mkdir -p "${DOWNLOAD_DIR}"
awk -F"\t" -vmax_dls="${MAX_DOWNLOADS}" -vdl_dir="${DOWNLOAD_DIR}" -vflist="$5" '
  FNR==NR {
    phenotypes[$1] = $1
  }
  ($5 == "variants.tsv.bgz") {
    split($6, arr, " ")
    if(system("[ ! -f " dl_dir "/../" arr[4] " ]") == 0)
      system(arr[1] " " arr[2] " " arr[3] " " dl_dir "/../" arr[4])
    else
      print dl_dir "/../" arr[4] " exists!"
  }
  ($5 ~ /male\.tsv\.bgz$/) {
    split($1, field, "_")
    if((field[1] in phenotypes) && (field[2] != "raw")){
      split($6, arr, " ")
      print "Processing " dl_dir "/" arr[4]
      if(system("[ ! -s " dl_dir "/" arr[4] " ]") == 0)
        system(arr[1] " " arr[2] " " arr[3] " " dl_dir "/" arr[4])
      else
        print dl_dir "/" arr[4] " exists!"
      dls += 1
      if ($5 ~ /\.male\.tsv\.bgz$/)
        print dl_dir "/" arr[4] >> flist ".male.txt"
      else
        print dl_dir "/" arr[4] >> flist ".female.txt"
    } else if (($1 ~ /^[^0-9]/) && ($1 != "N/A") && ($3 == "N/A")) {
      split($6, arr, " ")
      print "Processing " dl_dir "/" arr[4]
      if(system("[ ! -s " dl_dir "/" arr[4] " ]") == 0)
        system(arr[1] " " arr[2] " " arr[3] " " dl_dir "/" arr[4])
      else
        print dl_dir "/" arr[4] " exists!"
      dls += 1
      if ($5 ~ /\.male\.tsv\.bgz$/)
        print dl_dir "/" arr[4] >> flist ".male.txt"
      else
        print dl_dir "/" arr[4] >> flist ".female.txt"
    }
    if(dls == max_dls)
      exit
  }
' "${SELECTED_FIELDS}" <(gzip -dc "${MANIFEST}")
