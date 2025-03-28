#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=50000
#SBATCH --time=8:00:00
set -e
##***************************##
## RNAseq data preprocessing ##
##***************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 06.2022
## Creation: 06.2022
## Preprocessing of raw RNAseq data; from fastq to BAM


######------ Set environment ------######
module load python/2.7_intel
source activate my_ngs

print_help_message(){
    echo "
Usage: bash RNAseq_preprocessing.sh -d <sample-dir> -s <sample-name> -g <genome-dir> -a <annotation-file> -c <chrom-info> -o <output-dir>

We suggest to give full locations of input data.

Software required:
FastQC (version 0.11.9)
STAR (version STAR 2.6.0c)
samtools (version  1.12)
rseqc (version  2.6.4)
ucsc-wigtobigwig (version 377)
HTSeq (version 0.11.3)

Required:
    -d        |    Sample directory containing raw fastq files.

    -s        |    Sample (Li_DFxDF_++_2).

    -g        |    Genome directory including genomic index.

    -a        |    Annotation file for read count (gtf format).

    -c        |    File for BigWig generation (.chrom.sizes).

    -o        |    Output directory.

Misc:
    -h        |    Display this help message.
        "
    }

# cmd parsing: flags
while getopts :d:s:g:a:c:o: flag; do
    case $flag in
        d) RAWdir=$OPTARG;;
        s) SAMPLE=$OPTARG;;
        g) GENOMEdir=$OPTARG;;
        a) ANNOTATION=$OPTARG;;
        c) CHROMinfo=$OPTARG;;
        o) OUTPUTdir=$OPTARG;;
    esac
done
# check if required parameters are empty -> if so invoke help
[ -z $RAWdir ] && print_help_message && exit 1
[ -z $SAMPLE ] && print_help_message && exit 1
[ -z $GENOMEdir ] && print_help_message && exit 1
[ -z $ANNOTATION ] && print_help_message && exit 1
[ -z $CHROMinfo ] && print_help_message && exit 1
[ -z $OUTPUTdir ] && print_help_message && exit 1
# check if flag parameters are set
if [ ! -s $RAWdir ]; then
	echo "
Error: Fastq folder ${RAWdir} can not be found."
	exit 0
fi
if [ ! -s $GENOMEdir ]; then
	echo "
Error: Genome directory ${GENOMEdir} can not be found."
	exit 0
fi
if [ ! -s $ANNOTATION ]; then
	echo "
Error: Annotation file ${ANNOTATION} can not be found."
	exit 0
fi
if [ ! -s $CHROMinfo ]; then
	echo "
Error: Fastq folder ${CHROMinfo} can not be found."
	exit 0
fi

# make log file, rm if present
log=$OUTPUTdir"/"$SAMPLE".log"
exec > >(tee -i "$log")
exec 2>&1

echo -e '\n|---- Preprocessing of RNAseq data ----|\n'
echo "Input parameters:"
echo -e 'Sample name:' $SAMPLE
echo -e 'Genome directory:' $GENOMEdir
echo -e 'Annotation file:' $ANNOTATION

# make directories
mkdir -p $OUTPUTdir"/"$SAMPLE
RESULTS=$OUTPUTdir"/"$SAMPLE
mkdir -p $RESULTS
mkdir -p $RESULTS"/STAR_alignment_SNPsplit"


######------ Alignment ------######
echo -e '\n|-- Read alignment using STAR --|'
echo -e '\nplease wait...'
# get input path/name of reads and technical replicates
fastq_F=`ls $RAWdir | grep $SAMPLE.*_1.fastq.gz | awk -v sd=$RAWdir -v ORS="," '{print sd"/"$1}' | sed 's/,$//'`
fastq_R=`ls $RAWdir | grep $SAMPLE.*_2.fastq.gz| awk -v sd=$RAWdir -v ORS="," '{print sd"/"$1}' | sed 's/,$//'`

# perform alignment with STAR
time STAR \
--genomeDir $GENOMEdir \
--readFilesIn $fastq_F $fastq_R --readFilesCommand gunzip -c --runThreadN 12 --outStd Log \
--alignIntronMax 100000 --outFilterIntronMotifs RemoveNoncanonical \
--outFileNamePrefix $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE"_"  --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate
echo -e '...done!'
resultFile=$RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE"_Aligned.sortedByCoord.out.bam"
mv $resultFile $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE".bam"

# index bam file
echo -e '\nCreating index bam file...'
samtools index $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE".bam"
echo -e '...done!'


######------ BigWig file ------######
echo -e '\nCreating BigWig file...'
bam2wig.py -i $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE".bam" -s ${CHROMinfo} -o $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE --skip-multi-hits
rm -f $RESULTS"/STAR_alignment_SNPsplit"/$SAMPLE".wig"
echo -e '...done!'
echo -e '\nYour RNAseq data' $SAMPLE 'is now ready for downstream analysis.'
