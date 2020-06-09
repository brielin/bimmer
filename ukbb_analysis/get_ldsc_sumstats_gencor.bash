#!/bin/bash
#
# This script downloads the separated male and female UK biobank summary
# statistics from the Neale lab dropbox by parsing the LDSC results
# and analysis manifest. It pulls all sumstats with Z-score above 4 and
# medium or high confidence.
#
# Args:
#   $1: Manifest location.
#   $2: LDSC results location.
#   $3: Desired download directory.
#   $4: Maximum number of files to download.
#   $5: File prefix to write male/female file lists to.


MANIFEST=$1
LDSC_RESULTS=$2
DOWNLOAD_DIR=$3
MAX_DOWNLOADS=$4

[ ! -d "${DOWNLOAD_DIR}" ] && mkdir -p "${DOWNLOAD_DIR}"
awk -F"\t" -vmax_dls="${MAX_DOWNLOADS}" -vdl_dir="${DOWNLOAD_DIR}" -vflist="$5" '
  FNR==NR {
    if(($9 == "z4" || $9 == "z7") && ($10 == "medium" || $10 == "high")){
      phenotypes[$1] = $1
      print $1
    }
    next
  }
  ($1 in phenotypes) && ($4 == "both_sexes")  {
    split($7, arr, " ")
    print "Processing " dl_dir "/" arr[4]
    if(system("[ ! -s " dl_dir "/" arr[4] " ]") == 0)
      system(arr[1] " " arr[2] " " arr[3] " " dl_dir "/" arr[4])
    else
      print dl_dir "/" arr[4] " exists!"
    dls += 1

    if(dls == max_dls)
      exit
  }
' <(gzip -dc "${LDSC_RESULTS}") <(gzip -dc "${MANIFEST}")
