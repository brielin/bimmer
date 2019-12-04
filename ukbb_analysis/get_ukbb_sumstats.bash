#!/bin/bash
#
# This script downloads the separated male and female UK biobank summary
# statistics from the Neale lab dropbox by parsing the LDSC results
# and analysis manifest. It pulls all sumstats with Z-score above 4 and
# medium or high confidence.

MANIFEST="UKBB_GWAS_Manifest_201807.tsv.gz"
LDSC_RESULTS="ukb31063_h2_topline.02Oct2019.tsv.gz"
MAX_DOWNLOADS=4
DOWNLOAD_DIR="data"

awk -F"\t" -vmax_dls="${MAX_DOWNLOADS}" '
  FNR==NR {
    if(($9 == "z4" || $9 == "z7") && ($10 == "medium" || $10 == "high")){
      phenotypes[$1] = $1
    }
    next
  }
  {
    if(($1 in phenotypes) && ($4 == "male" || $4 == "female")){
      print $6
      system($6)
      dls += 1
    }
    if(dls == max_dls){
      exit
    }
  }
' <(gzip -dc "${LDSC_RESULTS}") <(gzip -dc "${MANIFEST}")

[ ! -d "${DOWNLOAD_DIR}" ] && mkdir -p "${DOWNLOAD_DIR}"
mv *.{male,female}.tsv.bgz "${DOWNLOAD_DIR}"
