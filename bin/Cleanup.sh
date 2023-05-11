#!/bin/bash

###  Inputs  ###
wkdir=$1
outdir=$2
allele=$3

cd ${wkdir}

###  Code  ###

rm ${outdir}/${allele}*.vcf.gz*
rm ${outdir}/${allele}_Genotype.tsv

touch ${outdir}/${allele}.done

exit 0
