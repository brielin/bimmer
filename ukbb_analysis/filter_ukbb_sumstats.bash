#!/bin/bash
#
# This preprocesses the UKBB sumstats by filtering down to only SNPs that
# will be used in the network MR analysis. That is, SNPs with pvalue < val
# for some phenotype in the analysis directory.
#
# Args:
#   $1: Directory with sumstats files to search.

# This needs to be run to get the SNPS_TO_USE arugment first!
# cat "${clump_dir}/"*"/snps_to_use.txt" | sort | uniq > "${filtered_dir}/all.snps_to_use.txt"


SUMSTATS_FILE=$1
VARIANTS_FILE=$2
SNPS_TO_USE=$3

raw_dir=$(dirname "${SUMSTATS_FILE}")
file_name=$(basename "${SUMSTATS_FILE}")
pheno_name="${file_name%%.tsv.bgz}"
clump_dir=$(dirname "${raw_dir}")"/clumped"
filtered_dir=$(dirname "${raw_dir}")"/filtered"

if [ ! -d "${filtered_dir}" ]; then
    mkdir "${filtered_dir}"
fi

echo "Processing ${pheno_name}."
if [ -f "${filtered_dir}/${pheno_name}_filtered.tsv.gz" ]; then
    echo "${filtered_dir}/${pheno_name}_filtered.tsv.gz already exists!"
else
    awk -F"\t" -vOFS="\t" '
      NR == FNR {
        snp_list[$1] = $1
        next
      }
      FNR == 1 {
        for(i=1; i<=NF; ++i){
          fields[$i] = i
        }
        print "rsid", "chr", "pos", "tstat", "n_complete_samples", "pval", "low_confidence_variant"
        next
      }            
      $(fields["rsid"]) in snp_list {
        print $(fields["rsid"]), $(fields["chr"]), $(fields["pos"]), $(fields["tstat"]), $(fields["n_complete_samples"]), $(fields["pval"]), $(fields["low_confidence_variant"])
      }
    ' "${SNPS_TO_USE}" <(paste <(zcat "${VARIANTS_FILE}") <(zcat "${SUMSTATS_FILE}" | cut -f 4-)) | gzip -c > "${filtered_dir}/${pheno_name}_filtered.tsv.gz"
fi
