##******************##
## Allelome.PRO run ##
##******************##
## Project: DissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 05.2022
## Submit samples for Allelome.PRO run to LRZ-cluster


######------ Set environment ------######
pipeline="./Allelome.PRO.sh"
input="./"
snps="./dissect_X_SNPs.bed"
output_dir="./05_AP_pseudobulk/"


######------ Set annotation ------######
# gene annotation
annotation="./annotation_us.bed"


######------ Set samples ------######
samples=("KO_bl6_Xa" "KO_cast_Xa" "WT_bl6_Xa" "WT_cast_Xa")


######------ Submit to LRZ ------######
for sample in "${samples[@]}"
do
  echo "Processing:"
  sbatch $pipeline -i $input${sample}".bam" -a $annotation -s $snps -r 1 -t 1 -o $output_dir
done
