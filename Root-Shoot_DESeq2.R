########################
#         ROOT        #
########################

setwd("~/RScounts/root/")

#sampleTable and factors
sampleTable_root=read.csv("~/RScounts/root/root_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_root$Genotype = factor(sampleTable_root$Genotype)
sampleTable_root$Gravity = factor(sampleTable_root$Gravity)
sampleTable_root$Run = factor(sampleTable_root$Run)
sampleTable_root$Position = factor(sampleTable_root$Position)
sampleTable_root$CVG=factor(paste(sampleTable_root$Gravity, sampleTable_root$Genotype, sep="."))

#remove GC samples
sampleTable_root=sampleTable_root[-c(21:24),]

#list of files
files.edgeR=as.character(sampleTable_root$Filename)

#create Edge object (remove low reads; goes down from 37334 genes to 19921)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
edge_root_counts=EDGE_filtered$counts

#normfactors
EDGE_N=calcNormFactors(EDGE_filtered)

#replace dispersions with model
Run=factor(sampleTable_root$Run)
Genotype=factor(sampleTable_root$Genotype)
Gravity=factor(sampleTable_root$Gravity)
CVG=factor(sampleTable_root$CVG)

design = model.matrix(~Run + Genotype + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)

#perform quasi-likelihood F-test
#(better with replicates, else use LRT)

fitQLF = glmQLFit(disp, design, robust=TRUE)

my.contrasts=makeContrasts()

QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)

FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05)

write.csv(topQLF.df, file="~/RScounts/root/edgeR_root_1291.csv")



########################
#         SHOOT        #
########################

setwd("~/RScounts/shoot/")

#sampleTable and factors
sampleTable_shoot=read.csv("~/RScounts/shoot/shoot_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_shoot$Genotype = factor(sampleTable_shoot$Genotype)
sampleTable_shoot$Gravity = factor(sampleTable_shoot$Gravity)
sampleTable_shoot$Run = factor(sampleTable_shoot$Run)
sampleTable_shoot$Position = factor(sampleTable_shoot$Position)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))

#remove GC samples
sampleTable_shoot=sampleTable_shoot[-c(15:20),]

#list of files
files.edgeR=as.character(sampleTable_shoot$Filename)

#create Edge object (remove low reads; goes down from 37334 genes to 20020)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
edge_shoot_counts=EDGE_filtered$counts

#normfactors
EDGE_N=calcNormFactors(EDGE_filtered)

#replace dispersions with model
Run=factor(sampleTable_shoot$Run)
Genotype=factor(sampleTable_shoot$Genotype)
Gravity=factor(sampleTable_shoot$Gravity)
CVG=factor(sampleTable_shoot$CVG)

design = model.matrix(~Run + Genotype + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)

#perform quasi-likelihood F-test
#(better with replicates, else use LRT)

fitQLF = glmQLFit(disp, design, robust=TRUE)

my.contrasts=makeContrasts()

QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)

FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05)

write.csv(topQLF.df, file="~/RScounts/shoot/edgeR_shoot_2539.csv")


############################################

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
files.edgeR=as.character(sampleTable_root$Filename)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
EDGE_N=calcNormFactors(EDGE_filtered)
Run=factor(sampleTable_root$Run)
Gravity=factor(sampleTable_root$Gravity)
CVG=factor(sampleTable_root$CVG)
design = model.matrix(~Run + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)
fitQLF = glmQLFit(disp, design, robust=TRUE)
QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)
FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05) #543
write.csv(topQLF.df, file="~/RScounts/root/edgeR_root_WT_543.csv")

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
files.edgeR=as.character(sampleTable_root$Filename)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
EDGE_N=calcNormFactors(EDGE_filtered)
Run=factor(sampleTable_root$Run)
Gravity=factor(sampleTable_root$Gravity)
CVG=factor(sampleTable_root$CVG)
design = model.matrix(~Run + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)
fitQLF = glmQLFit(disp, design, robust=TRUE)
QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)
FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05) #53
write.csv(topQLF.df, file="~/RScounts/root/edgeR_root_TR_53.csv")


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
sampleTable_shoot$Position = factor(sampleTable_shoot$Position)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))
files.edgeR=as.character(sampleTable_shoot$Filename)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
EDGE_N=calcNormFactors(EDGE_filtered)
Run=factor(sampleTable_shoot$Run)
Gravity=factor(sampleTable_shoot$Gravity)
CVG=factor(sampleTable_shoot$CVG)
design = model.matrix(~Run + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)
fitQLF = glmQLFit(disp, design, robust=TRUE)
QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)
FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05) #789
write.csv(topQLF.df, file="~/RScounts/shoot/edgeR_shoot_WT_789.csv")

## # # # # # # # # # # # 
# # #    TR only   # # #
# # # # # # # # # # # ##

setwd("~/RScounts/shoot/")
sampleTable_shoot=read.csv("~/RScounts/shoot/shoot_sampleTable.csv", sep=",", header=TRUE, stringsAsFactors = F)
sampleTable_shoot=sampleTable_shoot[sampleTable_shoot$Genotype == "TR" & sampleTable_shoot$Gravity != "GC",]
sampleTable_shoot$Genotype = factor(sampleTable_shoot$Genotype)
sampleTable_shoot$Gravity = factor(sampleTable_shoot$Gravity)
sampleTable_shoot$Run = factor(sampleTable_shoot$Run)
sampleTable_shoot$Position = factor(sampleTable_shoot$Position)
sampleTable_shoot$CVG=factor(paste(sampleTable_shoot$Gravity, sampleTable_shoot$Genotype, sep="."))
files.edgeR=as.character(sampleTable_shoot$Filename)
EDGE=readDGE(files=files.edgeR, columns=c(1,2))
EDGE_filtered <- EDGE[rowSums(1e+06 * EDGE$counts/expandAsMatrix(EDGE$samples$lib.size, dim(EDGE)) > 1) >= 2, ]
EDGE_N=calcNormFactors(EDGE_filtered)
Run=factor(sampleTable_shoot$Run)
Gravity=factor(sampleTable_shoot$Gravity)
CVG=factor(sampleTable_shoot$CVG)
design = model.matrix(~Run + Gravity)
disp = estimateDisp(EDGE_N, design, robust=TRUE)
fitQLF = glmQLFit(disp, design, robust=TRUE)
QLF = glmQLFTest(fitQLF)
topQLF = topTags(QLF, sort.by="PValue", p.value=0.05, n=Inf)
topQLF.df=as.data.frame(topQLF)
FDR = p.adjust(QLF$table$PValue, method="BH")
sum(FDR<0.05) #442
write.csv(topQLF.df, file="~/RScounts/shoot/edgeR_shoot_TR_442.csv")
