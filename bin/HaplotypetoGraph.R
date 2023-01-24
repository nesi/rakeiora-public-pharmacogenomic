####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################
#module purge; module load R/3.6.0-foss-2019a; R

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

#args <- commandArgs()
print(args)
dir <- args[1]
outdir <- args[2]
allele <- args[3]

###   Parameters  ### 

#Oone <- yes
#print("#### Options ####")
#print(paste0("Option 1: ", Oone ))

###   Directory creation  ### 
#setwd(dir)
print(paste0("Directory Set to ", dir))

###      Packages         ### 

library(tidyverse,quietly=TRUE,warn.conflicts=FALSE)

print("Packages Loaded")

####################################################################################################################
###                                               Load in Data                                                  ####
####################################################################################################################

read.table(paste0(outdir,allele,"_Haplotype.tsv"), header = T,sep="\t")  -> Input

####################################################################################################################
###                                                Allele Freqs                                                 ####
####################################################################################################################

Input %>% group_by(Haplotype) %>% summarise(Individuals = n(), HaplotypeFrequency = sum(Genotype)/(n()*2)) -> Haplotype_Summarises

write.table(Haplotype_Summarises, file = paste0(outdir,allele,"_Haplotype_Summary.tsv"),sep="\t",row.names=F)

Input$Genotype %>% gsub(0,"Homozygous-Reference",.) %>% gsub(1,"Heterozygous",.) %>% gsub(2,"Homozygous-Alternative",.) %>% factor(.,levels=c("Homozygous-Reference","Heterozygous","Homozygous-Alternative")) -> Input$Genotype_Written

####################################################################################################################
###                                               Plot Haplotypes                                               ####
####################################################################################################################

for(Haplotype in unique(Input$Haplotype)){
    Input %>% filter(Haplotype == Haplotype) -> Input_Plot
    
    plot <- ggplot(data=Input_Plot,aes(x=Genotype_Written)) +
        scale_y_continuous(expand = c(0,0),limits=c(0,nrow(Input_Plot))) +
        geom_bar(stat="count") +
        ggtitle(Haplotype) +
        ylab("Genotype Count") + 
        xlab(paste0(Haplotype," Genotype")) +
        scale_x_discrete(drop=FALSE) +
        theme(
            panel.grid = element_blank(),
            plot.title = element_text(size=50, face="bold", hjust = 0.5, margin = margin(t = 40, r = 0, b = 50, l = 0)),
            axis.text.x = element_text(angle = 0,size=34,vjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.text.y = element_text(angle = 0,size=34, margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.title.y = element_text(size=38, margin = margin(t = 0, r = 30, b = 0, l = 20)),
            axis.title.x = element_text(size=38, margin = margin(t = 30, r = 0, b = 20, l = 0)),
            axis.ticks=element_blank(),
            plot.background = element_rect(fill = 'white', colour = 'white'),
            panel.border = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(fill="white",colour="white",linetype="solid"),
            axis.line = element_line(colour = "black")
        )
    ggsave(plot=plot,filename=paste0(outdir,allele,"_Haplotype_Distribution.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)
}    

####################################################################################################################
###                                                      Fin                                                    ####
####################################################################################################################

quit()
