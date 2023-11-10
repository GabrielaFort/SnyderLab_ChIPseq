### DiffBind
# GF 042723
# This script will take a sample sheet file called diffbind.csv, run differential
# analysis based on the "condition" column, output DEseq results, .bed files
# of significant enriched peaks in each condition, and important graphs


#Load necessary packages
library("DiffBind")
library(tidyverse)
library(BiocParallel)
library("profileplyr")

register(SerialParam())


args <- commandArgs(trailingOnly = TRUE)
sampleSheet = args[1]
mydir = args[2]
fdr = args[3]

#must have a sample file with this name in the working directory
sample<-read.csv(paste(mydir, sampleSheet, sep="/"))


#read in the sample sheet
dbobj<-dba(sampleSheet=sample)


dbobj1<-dba.count(dbobj, bUseSummarizeOverlaps=TRUE, bRemoveDuplicates=FALSE)


dbobj2<-dba.normalize(dbobj1)


dbobj3<-dba.contrast(dbobj2, categories=DBA_CONDITION, minMembers = 2)


dbobj4<-dba.analyze(dbobj3, method=DBA_ALL_METHODS)


dba.show(dbobj4, bContrasts = TRUE)

dba.report(dbobj4)

###Extract results from DEseq and EDGER
res_deseq <- dba.report(dbobj4, method=DBA_DESEQ2, contrast = 1, th=1)
res_edger <- dba.report(dbobj4, method=DBA_EDGER, contrast = 1, th=1)

# Write to file
out_deseq <- as.data.frame(res_deseq)
out_edger <- as.data.frame(res_edger)

write.table(out_edger, file="./diffbind_edger.txt", sep="\t", quote=F, row.names=F)
write.table(out_deseq, file="./diffbind_deseq2.txt", sep="\t", quote=F, row.names=F)

# Write bed files
# Create bed files for each condition keeping only significant peaks (p < 0.05)
# in this case, condition 1 has a negative fold change and condition 2 has a positive fold change
# where the condition you added first in your sheet is condition 1

cond2_deseq <- out_deseq %>% 
  filter(FDR < fdr & Fold > 0) %>% 
  select(seqnames, start, end)


cond2_edger <- out_edger %>%
  filter(FDR < fdr & Fold > 0) %>%
  select(seqnames, start, end)
	
# Write to file
write.table(cond2_deseq, file="condition2_enriched_deseq.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(cond2_edger, file="condition2_enriched_edger.bed", sep="\t", quote=F, row.names=F, col.names=F)

cond1_deseq <- out_deseq %>% 
  filter(FDR < fdr & Fold < 0) %>% 
  select(seqnames, start, end)

cond1_edger <- out_edger %>%
  filter(FDR < fdr & Fold < 0) %>%
  select(seqnames, start, end)


# Write to file
write.table(cond1_deseq, file="condition1_enriched_deseq.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(cond1_edger, file="condition1_enriched_edger.bed", sep="\t", quote=F, row.names=F, col.names=F)


### Now exporting a bunch of potentially useful graphs ##

#heatmap of all merged sites
pdf("heatmap_all.pdf")
dba.plotHeatmap(dbobj4)
dev.off()

#heatmap of only differentially bound sites
pdf("heatmap_diff.pdf")
dba.plotHeatmap(dbobj4, contrast=1)
dev.off()

#PCA plot using all merged sites
pdf("PCA_all.pdf")
dba.plotPCA(dbobj4,DBA_TREATMENT,label=DBA_ID)
dev.off()

#PCA plot using only differentially bound sites
pdf("PCA_diff.pdf")
dba.plotPCA(dbobj4,contrast=1,label=DBA_ID)
dev.off()

#Volcano plot
pdf("volcano.pdf")
dba.plotVolcano(dbobj4)
dev.off()

#Box Plot
pdf("boxplot.pdf")
dba.plotBox(dbobj4)
dev.off()



#Binding affinity heatmap using differentially bound sites
hmap<-colorRampPalette(c("red","black","green"))(n=13)

pdf("binding_heatmap.pdf")
dba.plotHeatmap(dbobj4,contrast=1,correlations=FALSE, scale="row",colScheme=hmap)
dev.off()




#Tornado plot of differential analysis marking gained or lost sites
profiles <- dba.plotProfile(dbobj4)

pdf("tornado_plot.pdf")
dba.plotProfile(profiles, matrices_color=circlize::colorRamp2(breaks=seq(0,30,length=11),colors=rev(brewer.pal(11,'RdYlBu'))))
dev.off()




#Tornado plot of differential analysis with replicates shown separately
profiles1<-dba.plotProfile(dbobj4, merge=NULL)

pdf("tornado_plot_reps.pdf")
dba.plotProfile(profiles1)
dev.off()

# Venn diagram of the two different analyses i.e. DEseq2 vs EdgeR
pdf("peakspertest.pdf")
dba.plotVenn(dbobj4,contrast=1,method=DBA_ALL_METHODS)
dev.off()
