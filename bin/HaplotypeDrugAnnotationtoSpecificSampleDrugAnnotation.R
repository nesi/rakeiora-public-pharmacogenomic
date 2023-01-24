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
drug <- args[4]
filenm <- args[5]

###   Parameters  ### 

#Oone <- yes
#print("#### Options ####")
#print(paste0("Option 1: ", Oone ))

###   Directory creation  ### 
setwd(dir)
print(paste0("Directory Set to ", dir))

###      Packages         ### 

if(!require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)){
    install.packages("tidyverse",quiet=TRUE)
    library(tidyverse,quietly=TRUE,warn.conflicts=FALSE)
}

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

### Filtered Annotations
HaplotypeAnnotations_DrugSplit %>% filter(Drug == drug) -> HaplotypeAnnotations_DrugSplit_Filter
Sample_Index_HapDrug %>% filter(Drug == drug) -> Sample_Index_HapDrug_Filter

left_join(Sample_Index_HapDrug_Filter,HaplotypeAnnotations_DrugSplit_Filter,by="Drug") -> Sample_Index_HapDrug_Annotations
write.table(Sample_Index_HapDrug_Annotations, file = paste0(outdir,allele,"_Haplotype_Annotation_Summary_Prescribed_",drug,".tsv"),sep="\t",row.names=F)

####################################################################################################################
###                                                    Plot                                                     ####
####################################################################################################################

Sample_Index %>% select(Sample) %>% pull() %>% unique() -> All_ID
Sample_Index_HapDrug %>% filter(Drug == drug) %>% select(Sample) %>% pull() -> Prescribed_ID
setdiff(All_ID,Prescribed_ID) -> Unprescribed_ID

Haplotypes$DrugPrescribed[Haplotypes$Sample %in% Prescribed_ID] <- paste0("Prescribed\n",str_to_title(drug))
Haplotypes$DrugPrescribed[Haplotypes$Sample %in% Unprescribed_ID] <- paste0("Not Prescribed\n",str_to_title(drug))

factor(Haplotypes$DrugPrescribed,levels=c(paste0("Prescribed\n",str_to_title(drug)),paste0("Not Prescribed\n",str_to_title(drug)))) -> Haplotypes$DrugPrescribed

Haplotypes$Genotype %>% gsub(0,"Homozygous-Reference",.) %>% gsub(1,"Heterozygous",.) %>% gsub(2,"Homozygous-Alternative",.) %>% factor(.,levels=c("Homozygous-Reference","Heterozygous","Homozygous-Alternative")) -> Haplotypes$Genotype_Written

plot <- ggplot(data=Haplotypes,aes(x=Genotype_Written)) +
    facet_wrap(~ DrugPrescribed,drop = FALSE) +
    scale_y_continuous(expand = c(0,0),limits=c(0,nrow(Haplotypes))) +
    geom_bar(stat="count") +
    ggtitle(paste0(allele," Haplotype Distribution by ", str_to_title(drug)," Prescription")) +
    ylab("Genotype Count") + 
    xlab(paste0(allele," Genotype")) +
    scale_x_discrete(drop=FALSE) +
    theme(
        panel.grid = element_blank(),
        plot.title = element_text(size=50, face="bold", hjust = 0.5, margin = margin(t = 40, r = 0, b = 50, l = 0)),
        axis.text.y = element_text(angle = 0,size=34, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 0,size=22,vjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=38, margin = margin(t = 0, r = 30, b = 0, l = 20)),
        axis.title.x = element_text(size=38, margin = margin(t = 30, r = 0, b = 20, l = 0)),
        axis.ticks=element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_text(size=38, margin = margin(t = 20, r = 0, b = 20, l = 0)),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill="white",colour="white",linetype="solid"),
        axis.line = element_line(colour = "black")
    )
ggsave(plot=plot,filename=paste0(outdir,allele,"_Haplotype_Distribution_",drug,".pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

####################################################################################################################
###                                                      Fin                                                    ####
####################################################################################################################

quit()
