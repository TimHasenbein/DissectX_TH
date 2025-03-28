#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH --time=5:00:00
set -e
##****************************##
## ATACseq data preprocessing ##
##****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 02.2022
## Creation: 02.2022
## Quality control of ATAC-Seq data


######------ Set environment ------######
module load python/2.7_intel
module load gatk/
module load perl
source activate my_ngs

# help message
print_help_message(){
    echo "
The script will use aligned, sorted and indexed ATACseq bam files as outputted by ATAC_alignment.sh
Steps performed:
## 1. Filter aligned reads:
      # Remove reads derived from mitochondrial genome
      # Remove mapping artifacts (fragments < 38 bp6 and fragments > 2 kb) + low quality reads (MAPQ < 20)
      # Remove blacklist regions
      # Mark duplicates
Output: A folder will be generated in the same directory as for the alignment pipeline

Usage: sbatch <path/script> -d <path/data> -s <name-sample> -b <path/blacklist>

We suggest to give full locations of input files.

Software required:
python (version 2.7)
bedtools (≥ version 2.20.1)
SAMtools (≥ version 0.1.19)
R (≥ version 3.1.0)
Perl (≥ version 5.20.0)
bedToBigBed (≥ version 377)
Deeptools (≥ version 3.3.0)
GATK (≥ version 3.8)


R packages required:
ATACseqQC
If not installed, the script will try to install them in the default R library path.

Required:
    -d        |    Input sample directory.

    -s        |    Sample name.

    -b        |    Blacklist location and name.

Misc:
    -h        |    Display this help message.

    "
  }

# cmd parsing: flags
while getopts :d:s:b: flag; do
    case $flag in
        d) dir=$OPTARG;;
        s) sample=$OPTARG;;
        b) blacklist=$OPTARG;;
    esac
done

# check if required parameters are empty -> if so invoke help
[ -z $dir ] && print_help_message && exit 1
[ -z $sample ] && print_help_message && exit 1
[ -z $blacklist ] && print_help_message && exit 1

# check if flag parameters are set
if [ ! -s $dir ]; then
	echo "
Error: Sample directory ${dir} can not be found."
	exit 0
fi
if [ ! -s $blacklist ]; then
	echo "
Error: Blacklist file ${blacklist} can not be found."
	exit 0
fi

# make folder
mkdir -p $dir$sample"/03_quality_control/"
mkdir -p $dir$sample"/04_peak_calling/"

# get variables
qc_dir=$dir$sample"/03_quality_control/"
out=$qc_dir$sample
bam_dir=$dir$sample$"/02_bowtie2_alignment/"
sample_path=$bam_dir$sample

# make log file, rm if present
log=$out"_QC.log"
exec > >(tee -i "$log")
exec 2>&1
echo -e '|---- Preprocessing ATACseq reads ----|\n'
echo -e 'Sample: '$sample


######------  Filter RNAseq reads ------######
echo -e '\nNumber of reads in raw bam file: '`samtools view -c $sample_path.bam`
# filtering BAM files to remove mitochondrial reads and other alignments of no interest
echo -e '\n|-- Removing reads from mt genome --|\n'
samtools view -h $sample_path.bam | grep -v chrM | samtools sort -O bam -o $out.rmChrM.bam
echo -e 'Number of reads in mt reads removed bam file: '`samtools view -c $out.rmChrM.bam`
echo -e '\n|-- Removing mapping artefacts (fragment length <=38bp|>=2000kb) and low quality reads (MAPQ <20) --|\n'
samtools view  -h $out.rmChrM.bam | awk  'BEGIN{FS=OFS="\t"} function abs(v) {return v < 0 ? -v : v} \
   /^@/ || abs($9) <= 2000 && abs($9) >= 38 && $5 >= 20 {print}' | \
   samtools view  -bu -h -  | \
   samtools sort  -l 5 -o $out.rmChrM.rmAF.bam  -O BAM  -
echo -e 'Number of reads after filtering bam file: '`samtools view -c $out.rmChrM.rmAF.bam`
# Indexing the BAM file
samtools index  $out.rmChrM.rmAF.bam


######------  Removing blacklist genes ------######
echo -e '\n|-- Removing blacklist genes --|\n'
bedtools intersect -abam $out.rmChrM.rmAF.bam -b $blacklist -v | samtools sort -o $out.rmChrM.rmAF.rmBl.bam -
samtools index $out.rmChrM.rmAF.rmBl.bam $out.rmChrM.rmAF.rmBl.bam.bai
echo -e 'Number of reads after filtering bam file: '`samtools view -c $out.rmChrM.rmAF.rmBl.bam`


######------  Removing duplicates ------######
echo -e '\n|-- Marking duplicates --|\n'
export _JAVA_OPTIONS=-Djava.io.tmpdir=$out
gatk MarkDuplicates \
        -I $out.rmChrM.rmAF.rmBl.bam  \
        -O $out.rmChrM.rmAF.rmBl.rmDup.bam  \
        --REMOVE_DUPLICATES TRUE \
        --CREATE_INDEX TRUE \
        -M $out.metrics
samtools sort -o $out.rmChrM.rmAF.rmBl.rmDup.bam  $out.rmChrM.rmAF.rmBl.rmDup.bam
samtools index  $out.rmChrM.rmAF.rmBl.rmDup.bam $out.rmChrM.rmAF.rmBl.rmDup.bam.bai
echo -e 'Number of reads after removing duplicates: '`samtools view -c $out.rmChrM.rmAF.rmBl.rmDup.bam`


######------  Generate bigWig file ------######
echo -e '\n|-- Generating bigwig file --|\n'
bam2wig.py -i $out.rmChrM.rmAF.rmBl.rmDup.bam -s "/dss/dssfs02/lwp-dss-0001/pr23fa/pr23fa-dss-0000/andergassen_lab/RNAseq/reference_files/STAR_genome_mm10/genome_dir_gencodeM25_mm10.201911/GRCm38.chrom.sizes" -o $out.rmChrM.rmAF.rmBl.rmDup.bam --skip-multi-hits
rm $out.rmChrM.rmAF.rmBl.rmDup.bam.wig


######------  Summary ------######
echo -e '\n|-- ATACseq quality control done: --|\n'
# final summary + remove files
echo -e 'Number of reads in raw bam file: '`samtools view -c $sample_path.bam`
echo -e 'Number of reads in mt reads removed bam file: '`samtools view -c $out.rmChrM.bam`
echo -e 'Number of reads after filtering bam file for artefacts and MAPQ: '`samtools view -c $out.rmChrM.rmAF.bam`
echo -e 'Number of reads after removing blacklist: '`samtools view -c $out.rmChrM.rmAF.rmBl.bam`
echo -e 'Number of reads after removing duplicates: '`samtools view -c $out.rmChrM.rmAF.rmBl.rmDup.bam`


######------ Peak calling ------######
echo -e '\n|---- Peak calling with MACS2 ----|\n'
macs2 callpeak  -t $out.rmChrM.rmAF.rmBl.rmDup.bam  -f BAMPE -n $sample  -g mm --broad --keep-dup all --outdir $dir$sample$"/04_peak_calling/"
echo -e 'done!\n'
