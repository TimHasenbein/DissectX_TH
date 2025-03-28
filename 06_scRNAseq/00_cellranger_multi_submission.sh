#!/bin/bash
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=56000
#SBATCH --time=48:00:00
set -e
##*****************************##
## scRNAseq data preprocessing ##
##*****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 02.2022
## Creation: 02.2022
## Preprocessing of raw scRNAseq data by using the cellranger multi pipeline


cellranger multi --id=dissectX_scLFD --csv=./config_dissectX.csv
