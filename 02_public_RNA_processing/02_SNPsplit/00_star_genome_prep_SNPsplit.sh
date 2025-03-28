#!/bin/bash
#SBATCH -J STAR_indexgenome
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH --time=5:00:00
#SBATCH --mail-user=tim.hasenbein@tum.de
set -e
module load python/2.7_intel
source activate my_ngs
##*********************************##
## STAR: Genome index for SNPsplit ##
##*********************************##
## Project:  dissectX
## Tim Hasenbein
## Last modification 03.2022
## Creation: 03.2022
## Creating a genome index for STAR alignment with N-masked genome
## Files required:
## N-Masked reference genome (fasta)
## Annotation file (GTF)

STAR --runThreadN 6 \
     --runMode genomeGenerate \
     --genomeDir ./genome_preparation_STAR \
     --genomeFastaFiles ./genome.N-masked.fa \
     --sjdbGTFfile ./gencode.vM25.primary_assembly.annotation.gtf  \
     --sjdbOverhang 100
