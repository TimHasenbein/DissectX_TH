##****************************##
## Strand rule identification ##
##****************************##
## Project: Allelome.Link
## Tim Hasenbein
## Last modification 12.2022
## Creation: 10.2021
## Used to “guess” how RNA-seq sequencing were configured, particulary how reads were stranded for strand-specific RNA-seq data, through comparing the “strandness of reads” with the “standness of transcripts”.
## The “strandness of reads” is determiend from alignment, and the “standness of transcripts” is determined from annotation.


######------ Run infer_experiment.py ------######
infer_experiment.py -i $BAM -r $BED


######------ Output ------######
#This is PairEnd Data
#Fraction of reads failed to determine: 0.0520
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.0134
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.9346
