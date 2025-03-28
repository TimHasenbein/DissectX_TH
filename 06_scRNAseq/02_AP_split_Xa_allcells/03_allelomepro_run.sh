#!/bin/sh
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_large
#SBATCH --qos=cm2_large
#SBATCH --nodes=64
#SBATCH --tasks-per-node=1
#SBATCH --time=2-00:00:00
set -e
##**********************##
## Allelome.PRO jobfarm ##
##**********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 04.2022
## Run Allelome.PRO for each cell individually on chromosome level


######------ Do for each sample folder ------######
######------ get file names  ------######
# first get names of all files to include in the jobfarm run
 ls > files.txt

######------ make task-list ------######
# generate the task-db file
pipeline="./src/Allelome.PRO.sh"
input="./"
annotation="./chr_annotation_mm10.bed"
snps="./dissect_X_SNPs.bed"
output_dir="./"
min_read=1
total_read=10
cat $input"files.txt" | while read line
do
  echo "$pipeline -i $input${line} -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
done > task-list.txt


######------ Merge the task files manually  ------######


######------ run this script with sbatch as jobfarm ------######
module load slurm_setup
module use -p /lrz/sys/share/modules/extfiles/jobfarming
module load jobfarm
# start jobfarm 
jobfarm start "./task-list.txt"
