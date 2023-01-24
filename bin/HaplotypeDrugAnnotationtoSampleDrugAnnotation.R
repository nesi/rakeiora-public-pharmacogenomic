####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################
#module purge; module load R/3.6.0-foss-2019a; R

rm(list=ls())

#args <- commandArgs()
args = commandArgs(trailingOnly=TRUE)
print(args)
dir <- args[1]
outdir <- args[2]
allele <- args[3]
filenm <- args[4]

###   Parameters  ### 

#Oone <- yes
#print("#### Options ####")
#print(paste0("Option 1: ", Oone ))

###   Directory creation  ### 
setwd(dir)
print(paste0("Directory Set to ", dir))

###      Packages         ### 

library(tidyverse,quietly=TRUE,warn.conflicts=FALSE)
#if(!require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)){
#    install.packages("tidyverse",quiet=TRUE)
#    library(tidyverse,quietly=TRUE,warn.conflicts=FALSE)
#}

print("Packages Loaded")

####################################################################################################################
###                                               Load in Data                                                  ####
####################################################################################################################

read.table(filenm, header = T,sep="\t")  -> Drugs
read.table(paste0(outdir,allele,"_Haplotype.tsv"), header = T,sep="\t")  -> Haplotypes
read.delim(paste0(outdir,allele,"_Haplotype_PharmGKBAnnotation.tsv"), header = T,sep="\t",quote="")  -> HaplotypeAnnotations
colnames(HaplotypeAnnotations) <- c("ID","StudyHaplotype","Gene","Drug","PMID","Phenotype","Significance","Notes","Sentence","Alleles","Population")

####################################################################################################################
###                                                Data Manipulation                                            ####
####################################################################################################################

left_join(Haplotypes,Drugs) -> Sample_Index
as.character(HaplotypeAnnotations$Drug) %>% gsub("\\\"","",. ) -> HaplotypeAnnotations$Drug
HaplotypeAnnotations %>% separate_rows(., Drug, convert = TRUE) %>% unique() -> HaplotypeAnnotations_DrugSplit
Sample_Index %>% filter(Drug %in% HaplotypeAnnotations_DrugSplit$Drug) -> Sample_Index_HapDrug
Sample_Index_HapDrug %>% filter(Genotype != 0) -> Sample_Index_HapDrug_NonRef
left_join(Sample_Index_HapDrug_NonRef,HaplotypeAnnotations_DrugSplit,by="Drug") -> Sample_Index_HapDrug_NonRef_Annotations

write.table(Sample_Index_HapDrug_NonRef_Annotations, file = paste0(outdir,allele,"_SampleHaplotype_Annotation_Summary.tsv"),sep="\t",row.names=F)

####################################################################################################################
###                                                      Fin                                                    ####
####################################################################################################################

quit()
