####################     Variables       ###########################

configfile: "config.yaml"

ALLELES          = config["ALLELES"]
DRUGS            = config["DRUGS"]

outdir           = config["outdir"]
logdir           = config["logdir"]
singdir          = config["singdir"]
resourcedir      = config["resourcedir"]

RSingularity     = config["RSingularity"]
BCFSingularity   = config["BCFSingularity"]

pharmgkb         = config["pharmgkb"]
pharmvar         = config["pharmvar"]
inputvcffile     = config["inputvcffile"]
druglistfile     = config["druglistfile"]
idlistfile       = config["idlistfile"]

####################     Initialise       ###########################

rule all:
    input:
        expand(outdir + '{allele}_{drug}_Report.html', allele= ALLELES, drug= DRUGS)
        ,expand(outdir + '{allele}.done', allele= ALLELES)


#####################     Module 2 - Genotype Manipulation       ###########################

rule makeAlleleFiles:
    input:
        vcf = resourcedir + inputvcffile
    output:
        o1 = outdir + "{allele}.vcf.gz"
    log: 
        logdir + "makeAlleleFiles.{allele}.log"
    shell:
        "cp {input.vcf} {output.o1} &> {log}"

# Split out haplotype variants from raw VCF based on pharmvar VCF,
# need to split a gene specific VCF out first -
# need to figure out how to manipulate ALLELES array into GENE array
rule GeneVCFtoGenotype:
    input:
        outdir + "{allele}.vcf.gz"
    output:
        outdir + "{allele}_Genotype.tsv"
    log: 
        logdir + "GeneVCFtoGenotype.{allele}.log"
    container:
        singdir + BCFSingularity
    shell:
        "bin/GeneVCFtoGenotype.sh $PWD {outdir} {pharmvar} {wildcards.allele} &> {log}"
     
rule GenotypetoHaplotype:
    input:
        outdir + "{allele}_Genotype.tsv"
    output:
        outdir + "{allele}_Haplotype.tsv"
    log: 
        logdir + "GenotypetoHaplotype.{allele}.log"
    shell:
        "bin/GenotypetoHaplotype.sh $PWD {outdir} {wildcards.allele} &> {log}"
        
rule HaplotypetoGraph:
    input:
        outdir + "{allele}_Haplotype.tsv"
    output:
        outdir + "{allele}_Haplotype_Distribution.pdf"
    log: 
        logdir + "HaplotypetoGraph.{allele}.log"
    container:
        singdir + RSingularity
    shell:
        "Rscript bin/HaplotypetoGraph.R $PWD {outdir} {wildcards.allele} &> {log}"

#####################     Module 3 - Haplotype by Drug Annotations     ###########################

rule makeIDListFile:
    input:
        vcf = resourcedir + inputvcffile
    output:
        o1 = outdir + idlistfile
    container:
        singdir + BCFSingularity
    shell:
        "bcftools query -l {input.vcf} > {output.o1}"

rule GetDrugList:
    input:
        i1 = pharmgkb
        ,i2 = outdir + idlistfile
    output:
        o1 = outdir + druglistfile
    log: 
        logdir + "GetDrugList.log"
    container:
        singdir + RSingularity
    shell:
        "Rscript bin/MockDrugData.R $PWD {input.i1} {input.i2} {output.o1} &> {log}"

rule HaplotypeDrugAnnotation:
    input:
        outdir + "{allele}_Haplotype.tsv"
        ,pharmgkb
    output:
        outdir + "{allele}_Haplotype_PharmGKBAnnotation.tsv"
    log: 
        logdir + "HaplotypeDrugAnnotation.{allele}.log"
    shell:
        "bin/HaplotypeDrugAnnotation.sh $PWD {outdir} {pharmgkb} {wildcards.allele} &> {log}"

rule HaplotypeDrugAnnotationtoSampleDrugAnnotation:
    input:
        outdir + "{allele}_Haplotype_PharmGKBAnnotation.tsv"
        ,i2 = outdir + druglistfile
    output:
        outdir + "{allele}_SampleHaplotype_Annotation_Summary.tsv"
    log: 
        logdir + "HaplotypeDrugAnnotationtoSampleDrugAnnotation.{allele}.log"
    container:
        singdir + RSingularity
    shell:
        "Rscript bin/HaplotypeDrugAnnotationtoSampleDrugAnnotation.R $PWD {outdir} {wildcards.allele} {input.i2} &> {log}"

rule HaplotypeDrugAnnotationtoSpecificDrugAnnotation:
    input:
        i1 = outdir + druglistfile
        ,i2 = outdir + "{allele}_Haplotype_PharmGKBAnnotation.tsv"
    output:
        outdir + "{allele}_Haplotype_Distribution_{drug}.pdf",
        outdir + "{allele}_Haplotype_Annotation_Summary_Prescribed_{drug}.tsv"
    log: 
        logdir + "HaplotypeDrugAnnotationtoSpecificDrugAnnotation.{allele}.{drug}.log"
    container:
        singdir + RSingularity
    shell:
        "Rscript bin/HaplotypeDrugAnnotationtoSpecificSampleDrugAnnotation.R $PWD {outdir} {wildcards.allele} {wildcards.drug} {input.i1} &> {log}"

#####################     Report and Cleanup      ###########################

rule report:
    input:
        outdir + "{allele}_Haplotype_Distribution_{drug}.pdf",
        outdir + "{allele}_Haplotype_Annotation_Summary_Prescribed_{drug}.tsv",
        outdir + "{allele}_Haplotype_Distribution.pdf"
    output:
        outdir + "{allele}_{drug}_Report.html"
    log: 
        logdir + "FinalReport.{allele}.{drug}.log"
    run:
        from snakemake.utils import report

        report("""
        An example report
        ===================================
        
        Stuff happened with {wildcards.allele} and {wildcards.drug}
        
        """, output[0], File=input[0] )

rule cleanup:
    input:
        outdir + "{allele}_Haplotype_Distribution.pdf"
    output:
        outdir + "{allele}.done"   
    log: 
        logdir + "Cleanup.{allele}.log"
    shell:
        "bin/Cleanup.sh $PWD {outdir} {wildcards.allele} &> {log}"
