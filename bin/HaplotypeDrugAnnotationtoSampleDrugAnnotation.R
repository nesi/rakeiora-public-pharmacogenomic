####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################

rm(list=ls())

#args <- commandArgs()
args = commandArgs(trailingOnly=TRUE)
print(args)
dir <- args[1]
outdir <- args[2]
allele <- args[3]
drug <- args[4]
presfile <- args[5]
dispfile <- args[6]

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
###                                                Data Manipulation                                            ####
####################################################################################################################


#this assumes that the key drug information is the first word which may not hold true for all entries but is needed for practicality 

groups <- c("Prescribed","Dispensed")
prefixes <- c("pres","disp")

for(group in groups) { 
	prefix <- prefixes[which(groups==group)]

	### Format Participant Drug Data
	drugs <- read.table(get(paste0(prefix,"file")), header = FALSE, sep=",", colClasses=rep("character" ,2), stringsAsFactors = FALSE) %>% unique() %>% `colnames<-`(c("ID","Drug"))
	#drugs <- drugs %>% separate_rows(., Drug, convert = TRUE, sep = ";") %>% separate_rows(., Drug, convert = TRUE, sep = " with ") %>% unique() #separate out ; delimitated (prescribed) or ' with ' delimitated (dispensed) entries onto separate lines
	drugs$Drug <- drugs$Drug %>% trimws(which = "both") %>% tolower() #%>% gsub(" .*", "", .) #convert to lowercase and removal all but the first word #gsub("([A-Za-z]+).*", "\\1" ,.)
	drugs$ID <- toupper(drugs$ID) #genotype data has uppercase [A-Z], so converted here for consistency
	
	### Format Participant Haplotype Data
	Haplotypes <- read.table(paste0(outdir,allele,"_Haplotype.tsv"), header = T,sep="\t")
	Haplotypes$Genotype %>% gsub(0,"Homozygous-Reference",.) %>% gsub(1,"Heterozygous",.) %>% gsub(2,"Homozygous-Alternative",.) %>% factor(.,levels=c("Homozygous-Reference","Heterozygous","Homozygous-Alternative")) -> Haplotypes$Genotype_Written
	Haplotypes$ID <- toupper(Haplotypes$ID) #converted here for consistency
	
	### Merge Participant Drug and Haplotype Data
	Sample_Index <- left_join(Haplotypes,drugs)
	write.table(Sample_Index, file = paste0(outdir,allele,"_Haplotype_Annotation_Summary_",group,".tsv"),sep="\t",row.names=F)
	
	### Format PharmGKB Annotation Data
	HaplotypeAnnotations <- read.delim(paste0(outdir,allele,"_Haplotype_PharmGKBAnnotation.tsv"), header = T,sep="\t",quote="") %>% `colnames<-`(c("VariantAnnotationID","AnnotationHaplotype","Gene","Drug","PMID","Phenotype","Significance","Notes","Sentence","Alleles","Population"))
	HaplotypeAnnotations <- HaplotypeAnnotations %>% separate_rows(., Drug, convert = TRUE, sep = ",\"") %>% unique() #separate out , delimitated entries onto separate lines
	HaplotypeAnnotations$Drug <- HaplotypeAnnotations$Drug %>% trimws(which = "both") %>% as.character() %>% gsub("\\\"","",. ) %>% gsub(" .*", "", .) %>% tolower() #convert to lowercase, remove forward slahes, and removal all but the first word

	HaplotypeAnnotations %>% filter(Drug == drug) -> HaplotypeAnnotations_Filter
	write.table(HaplotypeAnnotations_Filter, file = paste0(outdir,allele,"_Haplotype_PharmGKBAnnotation_",drug,".tsv"),sep="\t",row.names=F)

	### Merged Participant Data Filtering
	Sample_Index %>% filter(Drug %in% HaplotypeAnnotations$Drug) -> Sample_Index_HapDrug
	Sample_Index_HapDrug %>% filter(Drug == drug) -> Sample_Index_HapDrug_Filter
	
	left_join(Sample_Index_HapDrug_Filter,HaplotypeAnnotations_Filter,by="Drug") -> Sample_Index_HapDrug_Annotations
	write.table(Sample_Index_HapDrug_Annotations, file = paste0(outdir,allele,"_Haplotype_Annotation_Summary_",group,"_",drug,".tsv"),sep="\t",row.names=F)
}

####################################################################################################################
###                                                    Plot                                                     ####
####################################################################################################################

for(group in groups) { 
	prefix <- prefixes[which(groups==group)]

	Sample_Index %>% select(ID) %>% pull() %>% unique() -> All_ID
	Sample_Index_HapDrug %>% filter(Drug == drug) %>% select(ID) %>% pull() -> Hit_ID
	setdiff(All_ID,Hit_ID) -> Miss_ID
	
	Haplotypes$DrugStatus <- NA
	Haplotypes$DrugStatus[Haplotypes$ID %in% Hit_ID] <- paste0(group,"\n",str_to_title(drug))
	Haplotypes$DrugStatus[Haplotypes$ID %in% Miss_ID] <- paste0("Not ",group,"\n",str_to_title(drug))
	
	factor(Haplotypes$DrugStatus,levels=c(paste0(group,"\n",str_to_title(drug)),paste0("Not ",group,"\n",str_to_title(drug)))) -> Haplotypes$DrugStatus
		
	plot <- ggplot(data=Haplotypes,aes(x=Genotype_Written)) +
	    facet_wrap(~ DrugStatus,drop = FALSE) +
	    scale_y_continuous(expand = c(0,0),limits=c(0,nrow(Haplotypes))) +
	    geom_bar(stat="count") +
	    ggtitle(paste0(allele," Haplotype Distribution by ", str_to_title(drug)," Status")) +
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
	ggsave(plot=plot,filename=paste0(outdir,allele,"_Haplotype_Distribution_",group,"_",drug,".pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)
}

####################################################################################################################
###                                                      Fin                                                    ####
####################################################################################################################

quit()
