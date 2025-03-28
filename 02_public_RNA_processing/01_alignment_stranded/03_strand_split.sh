#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=20000
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
sample_in_dir=$RAWdir$bam"/results/02_STAR_alignment/$bam"
outdir=$RAWdir$bam"/results/02_STAR_alignment/"


######------ Run strand separation script ------######
perl $sep_perl $sample_in_dir".bam" $strandrule $outdir


######------ Generate bigWigFiles ------######
# index if not indexed yet
if [ ! -s $sample_in_dir"_fwd.bam.bai" ]; then
  samtools sort $sample_in_dir"_fwd.bam" -o $sample_in_dir"_fwd_sort.bam"
	samtools index $sample_in_dir"_fwd_sort.bam"
fi
if [ ! -s $sample_in_dir"_rev.bam.bai" ]; then
  samtools sort $sample_in_dir"_rev.bam" -o $sample_in_dir"_rev_sort.bam"
	samtools index $sample_in_dir"_rev_sort.bam"
fi
# generate wig file
bam2wig.py --input-file=$sample_in_dir"_fwd_sort.bam"  --chromSize=$chrom_info --out-prefix=$sample_in_dir"_fwd_sort.bam"
bam2wig.py --input-file=$sample_in_dir"_rev_sort.bam"  --chromSize=$chrom_info --out-prefix=$sample_in_dir"_rev_sort.bam"
rm $output*.wig
"
echo "Sample .bam files "${sample}" successfully split"
