---
title: "heatmap"
output: html_document
---

setwd("~/IonTorrent/")
library(RColorBrewer)

##################
# Load/prep data #
##################

#fully setup countTable and sampleTable
counts=(read.table("dds_counts.csv", sep=",", header = TRUE, check.names = F))
rownames(counts)=counts[,1]
counts=counts[-1]

sampleTable=read.table("SampleTable.txt", sep="")
sampleTable=sampleTable[-(9:10),]
colnames(sampleTable)<-c("SampleName","FileName","Gravity", "Position", "Genotype", "Run")

#genotype/gravity vector (4 buckets)
CVG<-factor(paste(sampleTable$Gravity,".",sampleTable$Genotype,sep=""), levels=c("1G.WT", "MG.WT", "1G.TR", "MG.TR")) 
sampleTable$CVG=CVG

#genotype/gravity/run vector (8 buckets)
GGR<-factor(paste(sampleTable$Gravity,".",sampleTable$Genotype,".R",sampleTable$Run,sep=""), 
            levels=c("1G.WT.R1", "1G.WT.R2",
                     "MG.WT.R1", "MG.WT.R2",
                     "1G.TR.R1", "1G.TR.R2",
                     "MG.TR.R1", "MG.TR.R2"),
            ordered = TRUE)
sampleTable$GGR=GGR

#Per-genotype vectors
GenoWT<-factor(paste(sampleTable$Gravity,".",sampleTable$Genotype,".R",sampleTable$Run,sep=""), 
            levels=c("1G.WT.R1", "1G.WT.R2",
                     "MG.WT.R1", "MG.WT.R2"))

GenoTR<-factor(paste(sampleTable$Gravity,".",sampleTable$Genotype,".R",sampleTable$Run,sep=""), 
            levels=c("1G.TR.R1", "1G.TR.R2",
                     "MG.TR.R1", "MG.TR.R2"))

#factors for genotype and gravity
Geno=factor(sampleTable$Genotype, levels=c("WT","TR"))
Gravity=factor(sampleTable$Gravity, levels=c("1G", "MG"))

#filter out rows with sd=0
counts_filtered <- counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0),]

#make ordered counttable for buckets - order by vector
counts_orderedGGR=counts_filtered[,order(GGR)]

counts_orderedGeno=counts_filtered[,order(Geno)]

counts_orderedGrav=counts_filtered[,order(Gravity)]

##load DE lists
comb_grav_results1164 <- read_csv("~/deseq/Results/IonTorrent/comb_grav_results1164.csv")
genes1164=comb_grav_results1164[1]

WTgravity_results_05_547 <- read_csv("heatmap/WT_TR_lists/WTgravity_results_05_547.csv")
genesWT547=WTgravity_results_05_547[1]

TRgravity_results423 <- read_csv("heatmap/WT_TR_lists/TRgravity_results423.csv")
genesTR423=TRgravity_results423[1]

#subset counttable genes and select relevant samples
counts1164=counts_orderedGGR[which(rownames(counts_orderedGGR)%in% genes1164$X1),]

countsTR423=counts_orderedGeno[which(rownames(counts_orderedGeno)%in%genesTR423$X1),]
countsTR423_TR=countsTR423[,c(11:20)]

countsWT547=counts_orderedGeno[which(rownames(counts_orderedGeno)%in%genesWT547$X1),]
countsWT547_WT=countsWT547[,c(1:10)]

######################
# heatmap generating #
######################

library(gplots)
library(ape)
install.packages("ggdendro")
library("ggdendro")

setwd("~/IonTorrent/heatmap/WT_TR_lists/")

## setup colors
rc1 <- colorRampPalette(colors = c("purple", "white"), space="Lab")(8)   
rc2 <- colorRampPalette(colors = c("white", "green"), space="Lab")(8)

cols <- c(rc1, rc2)

#setup breaks 
pos_1 = seq(1, 0.5, length.out= 3)
pos_2 = seq(0.5, 0.001, length.out= 5)[-1]
pos_3 = seq(0.001, 0, length.out= 3)[-1]

neg_3 = seq(0, -0.001, length.out= 3) [-1]
neg_2 = seq(-0.001, -0.5, length.out= 5)[-1]
neg_1 = seq(-0.5, -1.0, length.out= 3)[-1]

breaks=c(pos_1, pos_2, pos_3, neg_3, neg_2, neg_1)

###################
# heatmap command #
###################

TR423_heatmap=pheatmap(countsTR423,
         scale = "row",
         color= cols,
         breaks = breaks,
         width=16,
         height = 9,
         cluster_cols = F,
         cluster_rows = TRUE,
         main = "423TR",
         annotation_names_row = TRUE,
         border_color = NA,
         na_col = "#ffffff",
         legend = T,
         show_rownames = F,
         keep.dendro = T)

##Extract clusters of similar genes

#pull dendrogram information
TR423_hclust=TR423_split$tree_row

#cut dendrogram at 8 groups
TR423_cuts=as.data.frame(sort(cutree(TR423_hclust, k=8)))
write.csv(x = TR423_cuts, file = "~/IonTorrent/heatmap/WT_TR_lists/TR423_cuts.csv")

#find number of leaves in each cut
count(TR423_cuts[,1] == 8)
TR423_cuts_sums=c(25,163,96,24,84,10,11,10)
write.csv(TR423_cuts_sums, file="~/IonTorrent/heatmap/WT_TR_lists/TR423_cuts_sums.csv")

#plot dendrogram and add height line at cut
plot(TR423_hclust, labels = F, main = "TR423_cluster_dend")
abline(h=6.49, col="red", lty=2, lwd=2)
