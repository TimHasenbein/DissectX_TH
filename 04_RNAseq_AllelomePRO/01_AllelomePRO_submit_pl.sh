#!/bin/bash
set -e
##***********************##
## Allelome.PRO placenta ##
##***********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 08.2021
## Creation: 08.2021
## Submit placenta samples for Allelome.PRO run to LRZ-cluster


######------ Submit to LRZ ------######
for sample in "${samples[@]}"
do
  echo "Processing:"
  sbatch $pipeline -I $input${sample}"/results/02_STAR_alignment/"${sample}".bam" -A $annotation -S $snps -O $output_dir
done
