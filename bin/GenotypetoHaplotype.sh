#!/bin/bash

###  Preamble  ###

#   This script takes variants from a VCF file that make up a specific *allele haplotype and condenses them from variant call genotype to *allele haplotype genotype for downstream use. 
#   The major assumption of this condensation is that all variants are required to be alternative for one or both alleles to achieve a haplotype genotype of heterozygous or homozygous reference, respectively (i.e. 0/1 + 1/1 + 1/1 = 0/1)
#   This assumption may not necessary be true for all star alleles (i.e. only one of the variants is required for the phenotype), but this is a limitation of the data we have currently.
#   Additionally, this does not consider the possibility that alleles are in trans (i.e. 1|0 + 0|1 = 0/1 despite being 0/0 by the previous definition). There is no guarantee that variants will be in phasing distance and we don't have parentals.
#   Technical assumptions: GT fields only from fifth column onwards, only ref and alt (i.e. 0 and 1, but not 2+), that genotypes will be seperate by / or |, and missing information (i.e ./.) is currently being handled as 0/0). 

###  Inputs  ###
wkdir=$1
outdir=$2
allele=$3

cd ${wkdir}

###  Code  ###
genotype=${outdir}/${allele}_Genotype.tsv

echo -e "ID\tHaplotype\tGenotype" > ${outdir}/${allele}_Haplotype.tsv
head -n1 ${genotype} | cut -f 3- | tr '\t' '\n' | while read gt
do
    id=${gt}
    geno=$(awk -F'\t' 'BEGIN {FS="\t";} NR == 1 { for(i=1; i<=NF; i++) { if($i=="'${gt}'") { col=i; } } } NR > 1 { print $col; }' ${genotype} | awk '{split($0,array,"/|\\|"); print array[1]+array[2]}' | awk 'min=="" || $1 < min {min=$1}; END{ print min}')
    echo -e "${id}\t${allele}\t${geno}" >> ${outdir}/${allele}_Haplotype.tsv
done

exit 0
