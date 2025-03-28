##****************************##
## Strand rule identification ##
##****************************##
## Project: General
## Tim Hasenbein
## Last modification 10.2021
## Creation: 10.2021
## Identify the strandrule for the sequencing data


######------ Make bed file from GTF ------######
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $ANNOTATION | gtf2bed - > gencode.vM25_gm35612_LXR.annotation.bed # awk adds placeholder were the field transcript id is missing
BED="./GENCODE_M25GRCm38.p6_201911/gencode.vM25_gm35612_LXR.annotation.bed"


######------ Run infer_experiment.py ------######
infer_experiment.py -i $BAM -r $BED


######------ Output ------######
##This is PairEnd Data
##Fraction of reads failed to determine: 0.0512
##Fraction of reads explained by "1++,1--,2+-,2-+": 0.0059
##Fraction of reads explained by "1+-,1-+,2++,2--": 0.9429
