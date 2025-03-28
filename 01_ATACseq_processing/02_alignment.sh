#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH -n 10
#SBATCH --time=6:00:00
set -e
##*******************##
## ATACseq alignment ##
##*******************##
## Project: dissectX
## Tim Hasenbein
## Last modification 03.2022
## Creation: 11.2021
## ATAC-seq pipeline for data pre-procession


######------ Set environment ------######
module load python/2.7_intel
source activate my_ngs

# cmd parsing: flags
while getopts :d:s:g:c:o: flag; do
    case $flag in
        d) RAWdir=$OPTARG;;
        s) SAMPLE=$OPTARG;;
        g) GENOMEdir=$OPTARG;;
        c) CHROMinfo=$OPTARG;;
        o) OUTPUTdir=$OPTARG;;
    esac
done

# make log file, rm if present
log=$OUTPUTdir"/"$SAMPLE".log"
exec > >(tee -i "$log")
exec 2>&1

# make directories
mkdir -p $OUTPUTdir"/"$SAMPLE
RESULTS=$OUTPUTdir"/"$SAMPLE
mkdir -p $RESULTS"/01_fastqc/"
mkdir -p $RESULTS"/02_bowtie2_alignment"
mkdir -p $RESULTS"/03_peak_calling"

echo -e '\n|---- Preprocessing of ATACseq data ----|\n'
echo "Input parameters:"
echo -e 'Sample name:' $SAMPLE
echo -e 'Genome index:' $GENOMEdir


######------ FastQC report ------#####
echo -e '\n|-- Quality check of raw data --|'
echo -e '\nGenerating fastqc reports...\n'
fastqc "$RAWdir"/"$SAMPLE"*.fastq.gz -o $RESULTS"/01_fastqc"
echo -e '...done!'


######------ Alignment ------######
echo -e '\n|-- Read alignment using Bowtie2 --|'
echo -e '\nplease wait...'
bowtie2  --very-sensitive -p 8 -x $GENOMEdir  -1 "$RAWdir"/"$SAMPLE"*1.fastq.gz -2 "$RAWdir"/"$SAMPLE"*2.fastq.gz  \
  |  samtools view  -u  -  \
  |  samtools sort  -o $RESULTS"/02_bowtie2_alignment/"$SAMPLE".bam"  -
# index bam file
echo -e '\nCreating index bam file...'
samtools index $RESULTS"/02_bowtie2_alignment"/$SAMPLE".bam"
echo -e '...done!'


######------ BigWig file ------######
echo -e '\nCreating BigWig file...'
bam2wig.py -i $RESULTS"/02_bowtie2_alignment/"$SAMPLE".bam" -s ${CHROMinfo} -o $RESULTS"/02_bowtie2_alignment"/$SAMPLE --skip-multi-hits
echo -e '...done!'
rm -f $RESULTS"/02_bowtie2_alignment"/$SAMPLE".wig.compress"
echo -e '\n...done'
