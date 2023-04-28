#!/bin/bash

###  Preamble  ###

#   This script filters a VCF (full raw call atm) using PharmVar VCFs (i.e. their definitions) for user specified star alleles. 
#   Filtered haplotype variant VCFs are then outputted as a Genotype TSV
#   At this stage, I have not put in a check to ensure that all variants from PharmVar allele VCF are present in the target VCF.
#   Therefore, it is possible that the variant that isn't call is not being included in haplotype calculations
#   Additionally, I have no check to see if the output VCF is entirely empty 

###  Inputs  ###
wkdir=$1
outdir=$2
pharmvar=$3
allele=$4

#wkdir=
#pharmvar=
#rawcall=

cd ${wkdir}

###  Filter gene VCF for variants part of inputed star allele  ###
gene=$(echo ${allele} | sed 's/\*.*//g')
varmatch=$(echo ${allele} | sed 's/\*/_/g')
varvcf=${pharmvar}/${gene}/GRCh38/${varmatch}.vcf

[[ -d "${pharmvar}" ]] || {
    echo "PharmVar directory ${pharmvar} not present"
    exit 1
}

[[ -d "${pharmvar}/${gene}" ]] || {
    echo "Gene directory for ${gene} not present in PharmVar data"
    exit 1
}

[[ -d "${pharmvar}/${gene}/GRCh38" ]] || {
    echo "GRCh38 directory for ${gene} not present in PharmVar data"
    exit 1
}

[[ -f "$varvcf" ]] || {
    echo "Allele ${allele} for gene ${gene} not present in PharmVar data"
    exit 1
}

# We can assume that all vcf files in the pharmvar have been pre-zipped and indexed
#if [[ ! -f ${varvcf}.gz.tbi ]]
#then
#    bgzip -c ${varvcf} > ${varvcf}.gz
#    tabix -p vcf ${varvcf}.gz
#fi

bcftools index --tbi ${outdir}/${allele}.vcf.gz

bcftools filter --set-GTs . -R ${varvcf}.gz -o ${outdir}/${allele}_present.vcf.gz -O z ${outdir}/${allele}.vcf.gz
bcftools index --tbi ${outdir}/${allele}_present.vcf.gz 

bcftools filter --set-GTs . -T ^${outdir}/${allele}_present.vcf.gz -o ${outdir}/${allele}_missing.vcf.gz -O z  ${varvcf}.gz
bcftools index --tbi ${outdir}/${allele}_missing.vcf.gz 

bcftools merge --missing-to-ref -o ${outdir}/${allele}_Genotype.vcf.gz -O z ${outdir}/${allele}_present.vcf.gz ${outdir}/${allele}_missing.vcf.gz 
bcftools index --tbi ${outdir}/${allele}_Genotype.vcf.gz 

###  Output star allele VCF as a Genotype TSV for input into the next stage  ###
sample=$(zcat ${outdir}/${allele}.vcf.gz | grep '^#CHROM' | cut -f 10-)
echo -e "CHROM\tPOS\t${sample}" > ${outdir}/${allele}_Genotype.tsv
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' ${outdir}/${allele}_Genotype.vcf.gz >> ${outdir}/${allele}_Genotype.tsv

exit 0
