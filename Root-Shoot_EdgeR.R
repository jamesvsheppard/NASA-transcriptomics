########################
#         ROOT        #
########################

#remove metadata (last 5 lines) from all countsfiles
for (entry in (list.files(path = "~/RScounts/root/", full.names = F))){
  hold=read.table(entry, sep="\t", header=FALSE )
  hold2=hold[-c(37336:37341),]
  write.table(hold2, file=entry, row.names=F, col.names=F, sep="\t", quote=F)}
setwd("~/RScounts/root/")
sampleTable_root=read.csv("~/RScounts/root/root_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)

#remove GC samples
sampleTable_root=sampleTable_root[-c(21:24),]

#setup factors and design formula
sampleTable_root$Genotype = factor(sampleTable_root$Genotype)
sampleTable_root$Gravity = factor(sampleTable_root$Gravity)
sampleTable_root$Run = factor(sampleTable_root$Run)
sampleTable_root$Position = factor(sampleTable_root$Position)
sampleTable_root$CVG=factor(paste(sampleTable_root$Gravity, sampleTable_root$Genotype, sep="."))

designCVG=formula(~Run + CVG)
designRun=formula(~Run + Genotype + Gravity)

#Setup dds object and remove rowSums < 1
ddsrootRGG=DESeqDataSetFromHTSeqCount(sampleTable_root, directory="~/RScounts/root/", design=designRGG) #37336
ddsrootRGG=ddsrootRGG[(rowSums(counts(ddsrootRGG)) > 1)] #27243


ddsrootRGG=estimateSizeFactors(ddsrootRGG)
ddsrootRGG=estimateDispersions(ddsrootRGG)
ddsrootRGG=DESeq(ddsrootRGG)

vstrootRGG=varianceStabilizingTransformation(ddsrootRGG, blind=FALSE)

#Grab normalized counts table
counts_deseq_root=counts(ddsrootRGG, normalized = TRUE)


#results extraction
results_rootRGG=results(ddsrootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_rootRGG$padj<0.05, na.rm=T) #1654
keep=results_rootRGG[which(results_rootRGG$padj<0.05),]
write.csv(keep, file="rootRGG_1654.csv")


plotPCA(vstrootRGG, intgroup=c("Gravity")) + ggtitle("Root") + geom_label(aes(label = name))

############

  
#######################################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################################################################

########################
#         SHOOT        #
########################

#remove tail 5 lines
for (entry in (list.files(path = "~/RScounts/shoot/", full.names = F))){
  hold=read.table(entry, sep="\t", header=FALSE )
  hold2=hold[-c(37336:37341),]
  write.table(hold2, file=entry, row.names=F, col.names=F, sep="\t", quote=F)}
setwd("~/RScounts/shoot/")
sampleTable_shoot=read.csv("~/RScounts/shoot/shoot_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)

#remove GC samples
sampleTable_shoot=sampleTable_shoot[-c(15:20),]

#setup factors
sampleTable_shoot$Genotype = factor(sampleTable_shoot$Genotype)
sampleTable_shoot$Gravity = factor(sampleTable_shoot$Gravity)
sampleTable_shoot$Run = factor(sampleTable_shoot$Run)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))
Gravity=sampleTable_shoot$Gravity

designCVG=formula(~Run + CVG)
designRGG=formula(~Run + Genotype + Gravity)

#dds
ddsshootRGG=DESeqDataSetFromHTSeqCount(sampleTable_shoot, directory="~/RScounts/shoot/", design=designRGG) #37336
ddsshootRGG=ddsshootRGG[(rowSums(counts(ddsshootRGG)) > 1)] #27243

vstshootRGG=varianceStabilizingTransformation(ddsshootRGG, blind=FALSE)

ddsshootRGG=estimateSizeFactors(ddsshootRGG)
ddsshootRGG=estimateDispersions(ddsshootRGG)

ddsshootRGG=DESeq(ddsshootRGG)
counts_deseq_shoot=counts(ddsshootRGG, normalized = TRUE)

#results extraction
results_shootRGG=results(ddsshootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_shootRGG$padj<0.05, na.rm=T) #3224
keep=results_shootRGG[which(results_shootRGG$padj<0.05),]
write.csv(keep, file="shootRGG_3224.csv")

#PCA plotting
plotPCA(vstshootRGG, intgroup=c("Gravity")) + ggtitle("Shoot") + geom_label(aes(label = name))


#######################################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################################################################

#                                       SPLIT OUT WT / TR                                             #

#######################################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################################################################

##############################
###########################
########################
#         ROOT        #
########################
###########################
##############################

## # # # # # # # # # # # 
# # #    WT only   # # #
# # # # # # # # # # # ##
setwd("~/RScounts/root/")
sampleTable_root=read.csv("~/RScounts/root/root_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_root=sampleTable_root[sampleTable_root$Genotype == "WT" & sampleTable_root$Gravity != "GC",]
sampleTable_root$Genotype = factor(sampleTable_root$Genotype)
sampleTable_root$Gravity = factor(sampleTable_root$Gravity)
sampleTable_root$Run = factor(sampleTable_root$Run)
sampleTable_root$Position = factor(sampleTable_root$Position)
sampleTable_root$CVG=factor(paste(sampleTable_root$Gravity, sampleTable_root$Genotype, sep="."))
designRGG=formula(~Run + Gravity)
ddsrootRGG=DESeqDataSetFromHTSeqCount(sampleTable_root, directory="~/RScounts/root/", design=designRGG) #37336
ddsrootRGG=ddsrootRGG[(rowSums(counts(ddsrootRGG)) > 1)] #27243
ddsrootRGG=estimateSizeFactors(ddsrootRGG)
ddsrootRGG=estimateDispersions(ddsrootRGG)
ddsrootRGG=DESeq(ddsrootRGG)

vstrootRGG=varianceStabilizingTransformation(ddsrootRGG, blind=FALSE)

#results extraction
results_rootRGG=results(ddsrootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_rootRGG$padj<0.05, na.rm=T) #1055
keep=results_rootRGG[which(results_rootRGG$padj<0.05),]
write.csv(keep, file="rootRG_WT_1055.csv")


## # # # # # # # # # # # 
# # #    TR only   # # #
# # # # # # # # # # # ##
setwd("~/RScounts/root/")
sampleTable_root=read.csv("~/RScounts/root/root_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_root=sampleTable_root[sampleTable_root$Genotype == "TR" & sampleTable_root$Gravity != "GC",]
sampleTable_root$Genotype = factor(sampleTable_root$Genotype)
sampleTable_root$Gravity = factor(sampleTable_root$Gravity)
sampleTable_root$Run = factor(sampleTable_root$Run)
sampleTable_root$Position = factor(sampleTable_root$Position)
sampleTable_root$CVG=factor(paste(sampleTable_root$Gravity, sampleTable_root$Genotype, sep="."))
designRGG=formula(~Run + Gravity)
ddsrootRGG=DESeqDataSetFromHTSeqCount(sampleTable_root, directory="~/RScounts/root/", design=designRGG) #37336
ddsrootRGG=ddsrootRGG[(rowSums(counts(ddsrootRGG)) > 1)] #27243
ddsrootRGG=estimateSizeFactors(ddsrootRGG)
ddsrootRGG=estimateDispersions(ddsrootRGG)
ddsrootRGG=DESeq(ddsrootRGG)

vstrootRGG=varianceStabilizingTransformation(ddsrootRGG, blind=FALSE)

#results extraction
results_rootRGG=results(ddsrootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_rootRGG$padj<0.05, na.rm=T) #1654
keep=results_rootRGG[which(results_rootRGG$padj<0.05),]
write.csv(keep, file="rootRG_TR_438.csv")


##############################
###########################
########################
#         SHOOT       #
########################
###########################
##############################

## # # # # # # # # # # # 
# # #    WT only   # # #
# # # # # # # # # # # ##

setwd("~/RScounts/shoot/")
sampleTable_shoot=read.csv("~/RScounts/shoot/shoot_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_shoot=sampleTable_shoot[sampleTable_shoot$Genotype == "WT" & sampleTable_shoot$Gravity != "GC",]
sampleTable_shoot$Genotype = factor(sampleTable_shoot$Genotype)
sampleTable_shoot$Gravity = factor(sampleTable_shoot$Gravity)
sampleTable_shoot$Run = factor(sampleTable_shoot$Run)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))
designRGG=formula(~Run + Gravity)
ddsshootRGG=DESeqDataSetFromHTSeqCount(sampleTable_shoot, directory="~/RScounts/shoot/", design=designRGG) #37336
ddsshootRGG=ddsshootRGG[(rowSums(counts(ddsshootRGG)) > 1)] #27243
ddsshootRGG=estimateSizeFactors(ddsshootRGG)
ddsshootRGG=estimateDispersions(ddsshootRGG)
ddsshootRGG=DESeq(ddsshootRGG)

vstshootRGG=varianceStabilizingTransformation(ddsshootRGG, blind=FALSE)

#results extraction
results_shootRGG=results(ddsshootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_shootRGG$padj<0.05, na.rm=T) #1819
keep=results_shootRGG[which(results_shootRGG$padj<0.05),]
write.csv(keep, file="shootRG_WT_1819.csv")

## # # # # # # # # # # # 
# # #    TR only   # # #
# # # # # # # # # # # ##

setwd("~/RScounts/shoot/")
sampleTable_shoot=read.csv("~/RScounts/shoot/shoot_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_shoot=sampleTable_shoot[sampleTable_shoot$Genotype == "TR" & sampleTable_shoot$Gravity != "GC",]
sampleTable_shoot$Genotype = factor(sampleTable_shoot$Genotype)
sampleTable_shoot$Gravity = factor(sampleTable_shoot$Gravity)
sampleTable_shoot$Run = factor(sampleTable_shoot$Run)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))
designRGG=formula(~Run + Gravity)
ddsshootRGG=DESeqDataSetFromHTSeqCount(sampleTable_shoot, directory="~/RScounts/shoot/", design=designRGG) #37336
ddsshootRGG=ddsshootRGG[(rowSums(counts(ddsshootRGG)) > 1)] #27243
ddsshootRGG=estimateSizeFactors(ddsshootRGG)
ddsshootRGG=estimateDispersions(ddsshootRGG)
ddsshootRGG=DESeq(ddsshootRGG)

vstshootRGG=varianceStabilizingTransformation(ddsshootRGG, blind=FALSE)

#results extraction
results_shootRGG=results(ddsshootRGG, contrast=c("Gravity", "1g", "MG"))
sum(results_shootRGG$padj<0.05, na.rm=T) #1203
keep=results_shootRGG[which(results_shootRGG$padj<0.05),]
write.csv(keep, file="shootRG_TR_1203.csv")