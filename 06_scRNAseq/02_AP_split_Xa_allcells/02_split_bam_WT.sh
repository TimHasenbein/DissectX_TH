#!/bin/bash
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000
#SBATCH -n 10
#SBATCH --time=08:00:00
set -e
##***********************##
## Split scRNA bam files ##
##***********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 04.2022
## Split the bam file for female/male KO/WT into individual per cell bam files


######------ Set environment ------######
module load python/


######------ generate log file ------######
log=$bam"_"$cellfile".log"
exec > >(tee -i "$log")
exec 2>&1


######------ split bam file ------######
echo -e "Cellfile: $cellfile"
echo -e "Bamfile: $bam"
sinto filterbarcodes -b $bam -c  $cellfile
echo -e "done!"

