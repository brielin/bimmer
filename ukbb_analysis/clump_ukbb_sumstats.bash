#!/bin/bash
#
# This script clumps the UKBB sumstats, one file at a time using a provided
# reference bed file. The reference bed should be split by chromosome. That
# is, if the first agrgument is ".../reference_panel_chr" it expects to find
# files ".../reference_panel_chr{1..22}.{bed, bim, fam}".
#
# Args:

# This needs to be run to get the SNP lists before we proceed!
# for i in $(seq 1 22); do
#   cut -f 2 "${REFERENCE_BED}${i}.bim" | sort | uniq -d > "${SUMSTATS_DIR}/ukbb_chr${i}_dups.txt"
# done


P_THRESH=0.0001
REFERENCE_BED=$1
DUPS_PATTERN=$2
SUMSTATS_FILE=$3
VARIANTS_FILE=$4

raw_dir=$(dirname "${SUMSTATS_FILE}")
file_name=$(basename "${SUMSTATS_FILE}")
pheno_name="${file_name%%.tsv.bgz}"
clump_dir=$(dirname "${raw_dir}")"/clumped"

echo "Processing ${pheno_name}."
if [ ! -d "${clump_dir}/${pheno_name}" ]; then
    mkdir "${clump_dir}/${pheno_name}"
fi

if [ -f "${clump_dir}/${pheno_name}/snps_to_use.txt" ]; then
    echo "${clump_dir}/${pheno_name}/snps_to_use.txt already exists!"
else
    # paste <(zcat "${VARIANTS_FILE}") <(zcat "${SUMSTATS_FILE}" | cut -f 4-) |
    # 	awk -F"\t" -vfname="${clump_dir}/${pheno_name}/" -vOFS="\t" '
    #       $1 == "variant" {
    #         for(i=1; i<=NF; ++i){
    #           fields[$i] = i
    #         }
    #         curr_chr = 0
    #         next
    #       }
    #       $(fields["chr"]) != curr_chr {
    #         curr_chr = $(fields["chr"])
    #         print "Splitting chromosome", curr_chr
    #         curr_fn = fname"chr"curr_chr".tsv"
    #         print "rsid", "chr", "pos", "tstat", "n_complete_samples", "pval", "low_confidence_variant" > curr_fn
    #       }
    #       {
    #         print $(fields["rsid"]), $(fields["chr"]), $(fields["pos"]), $(fields["tstat"]), $(fields["n_complete_samples"]), $(fields["pval"]), $(fields["low_confidence_variant"]) > curr_fn
    #       }'
    
    for i in $(seq 1 22); do
	if [ -s "${clump_dir}/${pheno_name}/chr${i}.log" ]; then
	    echo "${clump_dir}/${pheno_name}/chr${i}.log exists and is not empty!"
	else
	    plink --bfile "${REFERENCE_BED}${i}" \
		  --clump "${clump_dir}/${pheno_name}/chr${i}.tsv" \
		  --clump-r2 0.05 \
		  --clump-kb 500 \
		  --clump-snp-field "rsid" \
		  --clump-field "pval" \
		  --clump-p1 "${P_THRESH}" \
		  --exclude "${SUMSTATS_DIR}/${DUPS_PATTERN}${i}.txt" \
		  --out "${clump_dir}/${pheno_name}/chr${i}"
	fi
    done
    awk 'NF && $1 != "CHR" {print $3}' "${clump_dir}/${pheno_name}/chr"*".clumped" > "${clump_dir}/${pheno_name}/snps_to_use.txt"
fi
