#!/bin/bash
#
# This preprocesses the UKBB sumstats by filtering down to only SNPs that
# will be used in the network MR analysis. That is, SNPs with pvalue < val
# for some phenotype in the analysis directory.

P_THRESH=0.00001
MIN_MAF=0.001
SUMSTATS_DIR="data"
SNP_FILE="snps_to_use.txt"

# The sumstats format from Neale lab has optional columns for binary traits.
# Therefore we track column names rather than numbers.
find "${SUMSTATS_DIR}" -type f -name "*.tsv.bgz" |
    xargs gzip -dc |
    awk -F"\t" -vmin_maf="${MIN_MAF}" \
	-vp_thresh="${P_THRESH}" \
	-vsnp_file="${SNP_FILE}" '
      $1 == "variant" {
        ++file_num
	print "Processing file", file_num
        for(i=1; i<=NF; ++i){
	  fields[$i] = i
	}
	next
      }
      {
        if(($(fields["minor_AF"]) > min_maf) &&
            ($(fields["low_confidence_variant"]) == "false") &&
            ($(fields["pval"]) < p_thresh)){
	  snps_to_use[$(fields["variant"])] = $(fields["variant"])
	}
      }
      END {
        for(snp in snps_to_use){
	  print snp > snp_file
	}
      }
    '

sort -o "${SNP_FILE}" "${SNP_FILE}"

# I would much rather do this with pipes and all files at the same time
# so I don't have to re-read the SNP list. Or even better, do it in
# one awk command with two passes over each file just holding the snp list
# in memory.
for file in "${SUMSTATS_DIR}/"*.tsv.bgz; do
    prefix="${file%%.tsv.bgz}"
    awk -F"\t" '
      NR == FNR {
        snp_list[$1] = $1
	next
      }
      {
        if($1 in snp_list){
	  print $0
	}
      }
    ' "${SNP_FILE}" <(gzip -dc "${file}") |
	gzip -c > "${prefix}_filtered.tsv.gz"
done
