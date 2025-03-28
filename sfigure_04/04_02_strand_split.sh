#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=20000
#SBATCH --time=5:00:00
set -e

##***************************************##
## Divide strands of strand-specific seq ##
##***************************************##
## Project: ageing
## Tim Hasenbein
## Last modification 11.2021
## Creation: 10.2021
## Splitting bam file into 2 strand specific bam files
module load python/2.7_intel
source activate my_ngs


######------ BAM parsing ------######
while getopts :b:d:r:p: flag; do
    case $flag in
        b) bam=$OPTARG;;
        d) RAWdir=$OPTARG;;
        r) strandrule=$OPTARG;;
        p) sep_perl=$OPTARG;;
    esac
done


######------ Set environment ------######
chrom_info="./GRCm38.chrom.sizes"
sample_in_dir=$RAWdir$bam"/results/02_STAR_alignment/$bam"
outdir=$RAWdir$bam"/results/02_STAR_alignment/"


######------ Run strand separation script ------######
perl $sep_perl $sample_in_dir".bam" $strandrule $outdir


######------ Generate bigWigFiles ------######
# index if not indexed yet
if [ ! -s $sample_in_dir"_fwd.bam.bai" ]; then
	samtools index $sample_in_dir"_fwd.bam"
fi
if [ ! -s $sample_in_dir"_rev.bam.bai" ]; then
	samtools index $sample_in_dir"_rev.bam"
fi
# generate wig file
bam2wig.py --input-file=$sample_in_dir".bam" --chromSize=${chrom_info} --out-prefix=${sample_in_dir} --skip-multi-hits --strand=${strandrule}
# generate big wig file
wigToBigWig -clip $sample_in_dir".Forward.wig" $chrom_info $sample_in_dir"_fwd.bw"
wigToBigWig -clip $sample_in_dir".Reverse.wig" $chrom_info $sample_in_dir"_rev.bw"
# remove redundant files
rm -f $sample_in_dir".Forward.wig"
rm -f $sample_in_dir".Reverse.wig"
echo "Sample .bam files "${sample}" successfully split"
