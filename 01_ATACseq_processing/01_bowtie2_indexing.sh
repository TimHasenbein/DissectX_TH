#!/bin/bash
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=50000
#SBATCH --time=8:00:00
set -e
##*************************##
## Genome indexing bowtie2 ##
##*************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 10.2021
## Creation: 10.2021


######------ Set environment ------######
module load python/2.7_intel
source activate ngs


######------ Genome indexing ------######
bowtie2-build  ./GENCODE_M25GRCm38.p6_201911/GRCm38.primary_assembly.genome.fa  bowtie2_index
