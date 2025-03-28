#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH -n 10
#SBATCH --time=05:00:00
set -e
##***********************##
## Split scRNA bam files ##
##***********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 07.2023
## Creation: 04.2022
## Split the bam file faccording to Xa Xi status


######------ Set environment ------######
module load python/


######------ make cellfile bam file ------######
log=$bam"_"$cellfile".log"
exec > >(tee -i "$log")
exec 2>&1
echo -e "Cellfile: $cellfile"
echo -e "Bamfile: $bam"
sinto filterbarcodes -b $bam -c $cellfile
echo -e "done!"




