#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH --time=02:00:00
set -e
##**********************##
## Intersect ATAC peaks ##
##**********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 03.2022
## Creation: 03.2022
## intersect ATAC peaks per organ and replicate
## count ATAC-Seq peaks
## Here shown for one tissue sample


######------ Set environment ------######
module load python/2.7_intel
source activate my_ngs


######------ intersect MACS2 bed files per tissue ------######
intersectBed -a $fem_1 \
             -b $fem_2 \
             > $output"female_broad.bed"
intersectBed -a $male_1 \
             -b $male_2 \
             > $output"male_broad.bed"


## intersect and count ATAC peaks across genomic window
######------ Set environment ------######
module load python/2.7_intel
source activate my_ngs


######------ Submit to LRZ ------######
intersectBed -a $annotation"100kb_swindow.bed" \
             -b $bed"ce_female_broad.bed" \
             -wa -c > $output"ce_female_count_broad.bed"
intersectBed -a $annotation"100kb_swindow.bed" \
             -b $bed"ce_male_broad.bed" \
             -wa -c > $output"ce_male_count_broad.bed"