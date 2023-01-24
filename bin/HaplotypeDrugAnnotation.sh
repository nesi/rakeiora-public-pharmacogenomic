#!/bin/bash

###  Preamble  ###

#   This script pulls out PharmGKB annotations based on inputted allele and drug combinations

###  Inputs  ###
wkdir=$1
outdir=$2
pharmgkb=$3 
allele=$4

cd ${wkdir}

###  Create Directory  ###

###  Code  ###
# Add an escape '\' to allow appropriate partial matching with awk, else run into issues with the '*' in the allele name being interpreted as a wildcard 
searchallele=$(echo ${allele} | sed 's/.*\*/\\*/g')
gene=$(echo ${allele} | sed 's/\*.*//g')

# Create output file with the source file header 
head -n1 ${pharmgkb} > ${outdir}/${allele}_Haplotype_PharmGKBAnnotation.tsv
# Add annotation entries that partially match the user-specified allele, this will include entries with multiple different alleles other than the user-specified and reference allele. 
# This partial matches entries with a ' +', '/', and 'end of cell' after it, while removing inappropriate matches (i.e. *26 form *2) this 'should' be all possibilities. 
# Filtering based on the "Alleles" field rather than the "Variant/Haplotypes" field as the latter includes all haplotypes discussed in the paper, but they are no necessarily all included in the experimentation that is presented (an observation, not sure if this is always the case)
awk -F'\t' 'BEGIN {FS="\t";} $3 == "'${gene}'" && ( $10 ~ /'${searchallele}' +/ || $10 ~ /'${searchallele}'\// || $10 ~ /'${searchallele}'$/ ) && $7 == "yes" { print $0 }' ${pharmgkb} >> ${outdir}/${allele}_Haplotype_PharmGKBAnnotation.tsv
# Do we want to only include entries that are significant? Entires are either 'yes','no','not 'stated', only around 50% are 'yes' -> ONLY 'YES'
# I can switch the filtering to colname based for redundancy in case column order changes, but this in turn runs into colname change issues -> FINE AS IS

# Remove redundant columns??? 

exit
