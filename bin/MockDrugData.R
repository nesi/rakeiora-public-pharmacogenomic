####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################
#module purge; module load R/3.6.0-foss-2019a; R

rm(list=ls())

#args <- commandArgs()
args = commandArgs(trailingOnly=TRUE)

#print(args)
#allele <- args[6]
dir <- args[1]
pharmgkb <- args[2]
id <- args[3]
filenm <- args[4]

###   Parameters  ### 

#Oone <- yes
#print("#### Options ####")
#print(paste0("Option 1: ", Oone ))

###   Directory creation  ### 
setwd(dir)
#print(paste0("Directory Set to ", analysisDir))

###      Packages         ### 

library(tidyverse,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library(data.table,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

#if(!require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)){
#    install.packages("tidyverse",quiet=TRUE)
#    library(tidyverse,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
#}
#if(!require(reshape2,quietly=TRUE,warn.conflicts=FALSE)){
#    install.packages("reshape2")
#    library(reshape2,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
#}
#if(!require(data.table,quietly=TRUE,warn.conflicts=FALSE)){
#    install.packages("data.table")
#    library(data.table,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
#}

#print("Packages Loaded")

####################################################################################################################
###                                               Load in Data                                                  ####
####################################################################################################################

fread(pharmgkb, header = T,sep="\t",quote="")  -> Input
fread(id, header = F,sep="\t") -> IDs
colnames(IDs) <- "Sample"

####################################################################################################################
###                                                Filter Entries                                               ####
####################################################################################################################

#gene <- allele %>% gsub("\\*.*","",.)
gene <- "CYP2C19"
Input %>% filter(Gene == gene & Significance == "yes") %>% select(`Drug(s)`) %>% separate_rows(., `Drug(s)`, convert = TRUE) %>% unique()  -> Drugs

#Input %>% select(`Drug(s)`) %>% separate_rows(., `Drug(s)`, convert = TRUE,sep=",\"") -> Drugs
#Drugs$`Drug(s)` %>% gsub("\\","",.,fixed=T) %>% gsub("\"","",.,fixed=T) %>% gsub("\n","",.,fixed=T) -> Drugs$`Drug(s)`
#unique(Drugs$`Drug(s)`) -> List
#write.table(List, file="Drugs.tsv",row.names=F,col.names=F)

####################################################################################################################
###                                               Plot Haplotypes                                               ####
####################################################################################################################

AgeProbs <- data.frame(c("30-40","40-50","50-60","60-70","70-80"),c(0.15,0.20,0.25,0.30,0.10))
colnames(AgeProbs) <- c("Age","Prob")

sample(AgeProbs$Age, size = nrow(IDs),replace=T, prob = AgeProbs$Prob) -> IDs$Age
IDs$DrugNum <- 0
for(i in seq(1:length(IDs$Sample))){
    if(IDs$Age[i] == "30-40"){
        IDs$DrugNum[i] <- sample(0:2,size=1)
    }
    if(IDs$Age[i] == "40-50"){
        IDs$DrugNum[i] <- sample(0:3,size=1)
    }
    if(IDs$Age[i] == "50-60"){
        IDs$DrugNum[i] <- sample(1:4,size=1)
    }
    if(IDs$Age[i] == "60-70"){
        IDs$DrugNum[i] <- sample(2:6,size=1)
    }
    if(IDs$Age[i] == "70-80"){
        IDs$DrugNum[i] <- sample(3:7,size=1)
    }
}

rm(list=ls(pattern="^df_*"))
for(i in seq(1:length(IDs$Sample))){
	drugnum <- IDs$DrugNum[i]
	sample <- IDs$Sample[i]
	age <- IDs$Age[i]

	if(drugnum!=0){
		sample(Drugs$`Drug(s)`, size=drugnum, replace=TRUE) -> drugsample
		data.frame(sample,drugsample) -> df
    	assign(paste0("df_",sample),df)
     }
}

ls(pattern="^df_*") -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) %>% filter(drugsample != 'n')-> Combined_DrugSample

colnames(Combined_DrugSample) <- c("Sample","Drug")

write.table(Combined_DrugSample, file=filenm, quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                                      Fin                                                    ####
####################################################################################################################

quit()
