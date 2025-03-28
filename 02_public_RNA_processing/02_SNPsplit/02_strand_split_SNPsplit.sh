#!/bin/bash
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=56000
#SBATCH --time=5:00:00
set -e
##***************************************##
## Divide strands of strand-specific seq ##
##***************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 03.2022
## Creation: 10.2021
## Splitting bam file into 2 strand specific bam files
module load python/2.7_intel
source activate my_ngs


######------ BAM parsing ------######
while getopts :b:d: flag; do
    case $flag in
        b) bam=$OPTARG;;
        d) RAWdir=$OPTARG;;
    esac
done


######------ Set environment ------######
chrom_info="./GRCm38.chrom.sizes"
sep_perl="./03_separate_BAM_strand.pl"
strandrule="1+-,1-+,2++,2--"
sample_in_dir=$RAWdir"/01_alignment/"$bam"/STAR_alignment_SNPsplit/"$bam
outdir=$RAWdir"/02_strand_split/"$bam"/"
mkdir  -p $outdir


######------ Run strand separation script ------######
perl $sep_perl $sample_in_dir".bam" $strandrule $outdir


######------ Generate bigWigFiles ------######
samtools sort $outdir$bam"_fwd.bam" -o $outdir$bam"_fwd_sort.bam"
samtools index $outdir$bam"_fwd_sort.bam"
samtools sort $outdir$bam"_rev.bam" -o $outdir$bam"_rev_sort.bam"
samtools index $outdir$bam"_rev_sort.bam"
# generate big wig file
bam2wig.py --input-file=$outdir$bam"_fwd_sort.bam"  --chromSize=$chrom_info --out-prefix=$outdir$bam"_fwd_sort.bam"
bam2wig.py --input-file=$outdir$bam"_rev_sort.bam"  --chromSize=$chrom_info --out-prefix=$outdir$bam"_rev_sort.bam"
rm $outdir$bam*.wig
echo "Sample .bam files "${sample}" successfully split"
