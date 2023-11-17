### DiffBind
# GF 042723
# This script will take a sample sheet file called diffbind.csv, run differential
# analysis based on the "condition" column, output DEseq results, .bed files
# of significant enriched peaks in each condition, and important graphs


# Install necessary packages
package_list <- c("DiffBind","tidyverse","BiocParallel","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}

# Need older version of ggtree installed to be compatible with profileplyr (necessary for tornado plots)
if(!"profileplyr" %in% installed.packages()) {
  library(devtools)
  devtools::install_github("YuLab-SMU/ggtree")
  BiocManager::install("profileplyr")
}

# Load packages
library("DiffBind")
library("tidyverse")
library("BiocParallel")
library("profileplyr")

register(SerialParam())


args <- commandArgs(trailingOnly = TRUE)
sampleSheet = args[1]
mydir = args[2]
fdr = args[3]

sample<-read.csv(paste(mydir, sampleSheet, sep="/"))


# Perform diffbind steps according to tutorial
dbobj<-dba(sampleSheet=sample)

dbobj<-dba.count(dbobj, bUseSummarizeOverlaps=TRUE, bRemoveDuplicates=FALSE)

dbobj<-dba.normalize(dbobj)

dbobj<-dba.contrast(dbobj, categories=DBA_CONDITION, minMembers = 2)

dbobj<-dba.analyze(dbobj, method=DBA_ALL_METHODS)

# Gives little summary
dba.show(dbobj, bContrasts = TRUE)

dba.report(dbobj)

###Extract results from both methods: DEseq and EDGER
db_results <- dba.report(dbobj, method=DBA_ALL_METHODS, contrast = 1, th=1)

# Write to file - bed file format
out_results <- as.data.frame(db_results)

# Write without column names to be compatible with bed file formatting - but include in documentation what each column is...
write.table(out_results, file="./diffbind_results.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Write bed files
# Create bed files for each condition keeping only significant peaks (according to user defined fdr)
# in this case, condition 1 has a negative fold change and condition 2 has a positive fold change
# where the condition you added first in your sheet is condition 1

cond1_diff <- out_results %>%
  filter(FDR < fdr & Fold < 0)

cond2_diff <- out_results %>% 
  filter(FDR < fdr & Fold > 0)

# Write to file
write.table(cond1_diff, file="condition1_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(cond2_diff, file="condition2_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


### Now exporting a bunch of potentially useful graphs ##

#heatmap of all merged sites
pdf("heatmap_all.pdf")
dba.plotHeatmap(dbobj)
dev.off()

#heatmap of only differentially bound sites
pdf("heatmap_diff.pdf")
dba.plotHeatmap(dbobj, contrast=1)
dev.off()

#PCA plot using all merged sites
pdf("PCA_all.pdf")
dba.plotPCA(dbobj,DBA_TREATMENT,label=DBA_ID)
dev.off()

#PCA plot using only differentially bound sites
pdf("PCA_diff.pdf")
dba.plotPCA(dbobj,contrast=1,label=DBA_ID)
dev.off()

#Volcano plot
pdf("volcano.pdf")
dba.plotVolcano(dbobj)
dev.off()

#Box Plot
pdf("boxplot.pdf")
dba.plotBox(dbobj)
dev.off()


#Binding affinity heatmap using differentially bound sites
hmap<-colorRampPalette(c("red","black","green"))(n=13)

pdf("binding_heatmap.pdf")
dba.plotHeatmap(dbobj,contrast=1,correlations=FALSE, scale="row",colScheme=hmap)
dev.off()


#Tornado plot of differential analysis marking gained or lost sites
profiles <- dba.plotProfile(dbobj)

pdf("tornado_plot.pdf")
dba.plotProfile(profiles)
dev.off()


#Tornado plot of differential analysis with replicates shown separately
profiles1<-dba.plotProfile(dbobj, merge=NULL)

pdf("tornado_plot_reps.pdf")
dba.plotProfile(profiles1)
dev.off()

# Venn diagram of the two different analyses i.e. DEseq2 vs EdgeR
pdf("peakspertest.pdf")
dba.plotVenn(dbobj,contrast=1,method=DBA_ALL_METHODS)
dev.off()


