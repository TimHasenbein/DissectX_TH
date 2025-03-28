#!/bin/bash
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=56000
#SBATCH --time=2:00:00
set -e
##****************##
##  Run SNPsplit  ##
##****************##
## Tim Hasenbein
## Last modification 11.2023
## Creation: 06.2022
## Run SNPsplit to split the sequencing reads to the alleles
module load python/2.7_intel
source activate my_ngs


######------ Set environment ------######
chrom_info="./GRCm38.chrom.sizes"
ref_snp="./all_CAST_EiJ_SNPs_FVB_NJ_reference.based_on_GRCm38.txt"
input="./02_strand_split/"
output="./03_snp_split/"


######------ run SNP-split-----######
# CxF
perl ./aBrain_CxF_fwd.bam
perl ./aBrain_CxF_rev.bam
# FxC
perl ./aBrain_FxC_fwd.bam
perl ./aBrain_FxC_rev.bam


######------ sort/index allele-specific files-----######
# CxF fwd
samtools sort $input"aBrain_CxF/aBrain_CxF_fwd.genome1.bam" -o $output"aBrain_CxF_fwd.genome1_sort.bam"
samtools sort $input"aBrain_CxF/aBrain_CxF_fwd.genome2.bam" -o $output"aBrain_CxF_fwd.genome2_sort.bam"
samtools index $output"aBrain_CxF_fwd.genome1_sort.bam"
samtools index $output"aBrain_CxF_fwd.genome2_sort.bam"
bam2wig.py --input-file=$output"aBrain_CxF_fwd.genome1_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_CxF_fwd.genome1_sort.bam" --skip-multi-hits
bam2wig.py --input-file=$output"aBrain_CxF_fwd.genome2_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_CxF_fwd.genome2_sort.bam" --skip-multi-hits
# CxF rev
samtools sort $input"aBrain_CxF/aBrain_CxF_rev.genome1.bam" -o $output"aBrain_CxF_rev.genome1_sort.bam"
samtools sort $input"aBrain_CxF/aBrain_CxF_rev.genome2.bam" -o $output"aBrain_CxF_rev.genome2_sort.bam"
samtools index $output"aBrain_CxF_rev.genome1_sort.bam"
samtools index $output"aBrain_CxF_rev.genome2_sort.bam"
bam2wig.py --input-file=$output"aBrain_CxF_rev.genome1_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_CxF_rev.genome1_sort.bam" --skip-multi-hits
bam2wig.py --input-file=$output"aBrain_CxF_rev.genome2_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_CxF_rev.genome2_sort.bam" --skip-multi-hits
# FxC fwd
samtools sort $input"aBrain_FxC/aBrain_FxC_fwd.genome1.bam" -o $output"aBrain_FxC_fwd.genome1_sort.bam"
samtools sort $input"aBrain_FxC/aBrain_FxC_fwd.genome2.bam" -o $output"aBrain_FxC_fwd.genome2_sort.bam"
samtools index $output"aBrain_FxC_fwd.genome1_sort.bam"
samtools index $output"aBrain_FxC_fwd.genome2_sort.bam"
bam2wig.py --input-file=$output"aBrain_FxC_fwd.genome1_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_FxC_fwd.genome1_sort.bam" --skip-multi-hits
bam2wig.py --input-file=$output"aBrain_FxC_fwd.genome2_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_FxC_fwd.genome2_sort.bam" --skip-multi-hits
# FxC rev
samtools sort $input"aBrain_FxC/aBrain_FxC_rev.genome1.bam" -o $output"aBrain_FxC_rev.genome1_sort.bam"
samtools sort $input"aBrain_FxC/aBrain_FxC_rev.genome2.bam" -o $output"aBrain_FxC_rev.genome2_sort.bam"
samtools index $output"aBrain_FxC_rev.genome1_sort.bam"
samtools index $output"aBrain_FxC_rev.genome2_sort.bam"
bam2wig.py --input-file=$output"aBrain_FxC_rev.genome1_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_FxC_rev.genome1_sort.bam" --skip-multi-hits
bam2wig.py --input-file=$output"aBrain_FxC_rev.genome2_sort.bam" --chromSize=$chrom_info --out-prefix=$output"aBrain_FxC_rev.genome2_sort.bam" --skip-multi-hits
rm $output*.wig