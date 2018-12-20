setwd("~/IonTorrent/")

#########################
# load and prepare data #
#########################

#load sampletable
sampleTable<-read.table("SampleTable.txt", sep = "\t", header=F)
colnames(sampleTable)<-c("SampleName","FileName","Gravity", "Position", "Genotype", "Run")
#remove GC samples
sampleTable=sampleTable[sampleTable$Gravity != "GC",]

#setup factors
sampleTable$Gravity=factor(sampleTable$Gravity)
sampleTable$Run=factor(sampleTable$Run)
sampleTable$Genotype=factor(sampleTable$Genotype)
sampleTable$Position=factor(sampleTable$Position)
Gravity = sampleTable$Gravity 

#setup re-ordered combinedvectorsgravity factor
CombinedvectorsGravity<-factor(paste(sampleTable$Gravity,".",sampleTable$Genotype,sep=""))
CVG = factor(CombinedvectorsGravity, levels=c("1G.WT", "MG.WT", "1G.TR", "MG.TR"), ordered = TRUE)
sampleTable$CVG=CVG
sampleTable$CVG=factor(sampleTable$CVG)

#load full count table
countTable=read.table("dds_counts.csv", sep=",", header=TRUE)
rownames(countTable)=countTable$X
countTable[1]=NULL
#remove GC samples and mapping statistic rows
countTable=as.data.frame(countTable)
countTable=countTable[-c(9:10)]
countTable=countTable[-c(33603:33607),]

#re-order counts table by CVG
orderedCounts=countTable[,order(CVG)]

#make copy of plot.table to preserve original
bucket.table=plot.table

###########################
# scale data for plotting #
###########################

#avg of 1GWT against which to "normalize"; (first four columns)
avg.1GWT=as.matrix(rowMeans(bucket.table[c(1:4)]))

#divide every cell by averaged 1GWT values
bucket.scaled=bucket.table/avg.1GWT

#remove rows where 1GWT avg=0, causing NaN value
bucket.scaled <- na.omit(bucket.scaled)

#average together all samples in each "bucket"
avg.1GWT=as.matrix(rowMeans(bucket.scaled[c(1:4)]))
avg.MGWT=as.matrix(rowMeans(bucket.scaled[c(5:10)]))
avg.1GTR=as.matrix(rowMeans(bucket.scaled[c(11:14)]))
avg.MGTR=as.matrix(rowMeans(bucket.scaled[c(15:20)]))

#combine the four "bucket average" vectors into one final matrix
bucket.final=as.data.frame(cbind(avg.1GWT, avg.MGWT, avg.1GTR, avg.MGTR))
colnames(bucket.final)=c("1G.WT", "MG.WT", "1G.TR", "MG.TR")
colnames=colnames(bucket.final)
bucket.final=as.matrix(bucket.final)

setwd("~/IonTorrent/composite/")

#load 5 archetype-lists
type1l=read.table(file="type1.csv", sep=",", header=F, stringsAsFactors = FALSE)
type2l=read.table(file="type2.csv", sep=",", header=F, stringsAsFactors = FALSE)
type3l=read.table(file="type3.csv", sep=",", header=F, stringsAsFactors = FALSE)
type4l=read.table(file="type4.csv", sep=",", header=F, stringsAsFactors = FALSE)
type5l=read.table(file="type5.csv", sep=",", header=F, stringsAsFactors = FALSE)

#subset main bucket table for five "types"
type1=as.data.frame(bucket.final[which(rownames(bucket.final)%in% type1l$V1),]) #17
type2=as.data.frame(bucket.final[which(rownames(bucket.final)%in% type2l$V1),]) #13
type3=as.data.frame(bucket.final[which(rownames(bucket.final)%in% type3l$V1),]) #8
type4=as.data.frame(bucket.final[which(rownames(bucket.final)%in% type4l$V1),]) #5
type5=as.data.frame(bucket.final[which(rownames(bucket.final)%in% type5l$V1),]) #5

#########################################################
##################    PLOTTING     ######################
#########################################################

#arrange data into table suitable for ggplot, plot, and save (type1) 14
type1.plot.table=c()
type1.plot.table=as.data.frame(
  cbind(
    c(type1$`1G.WT`,type1$MG.WT,type1$`1G.TR`,type1$MG.TR), 
    rep(rownames(type1),4), 
    c(rep("1G.WT",14),rep("MG.WT",14),rep("1G.TR",14),rep("MG.TR",14))
  ), stringsAsFactors=FALSE)
colnames(type1.plot.table)=c("scaled_cts", "gene", "condition")
type1.plot.table$scaled_cts=as.numeric(type1.plot.table$scaled_cts)
type1.plot.table$condition=factor(type1.plot.table$condition, levels = c("1G.WT", "MG.WT", "1G.TR", "MG.TR"))
b = ggplot(type1.plot.table, aes(x=condition, y=scaled_cts)) +
  scale_y_continuous(name="Normalized Counts (Scaled)", limits=c(0.9, 1.55)) +
  geom_jitter(position=position_jitter(0.3), aes(shape=gene)) +
  scale_shape_manual(values=c(1:14)) +
  theme_minimal()
ggsave("type1_composite.png", width=8, height=5)

#arrange data into table suitable for ggplot, plot, and save (type2) 16
type2.plot.table=c()
type2.plot.table=as.data.frame(
  cbind(
    c(type2$`1G.WT`,type2$MG.WT,type2$`1G.TR`,type2$MG.TR), 
    rep(rownames(type2),4), 
    c(rep("1G.WT",16),rep("MG.WT",16),rep("1G.TR",16),rep("MG.TR",16))
  ), stringsAsFactors=FALSE)
colnames(type2.plot.table)=c("scaled_cts", "gene", "condition")
type2.plot.table$scaled_cts=as.numeric(type2.plot.table$scaled_cts)
type2.plot.table$condition=factor(type2.plot.table$condition, levels = c("1G.WT", "MG.WT", "1G.TR", "MG.TR"))
b = ggplot(type2.plot.table, aes(x=condition, y=scaled_cts, group=gene)) +
  scale_y_continuous(name="Normalized Counts (Scaled)", limits=c(0.8, 1.4)) +
  geom_jitter(position=position_jitter(0.3), aes(shape=gene)) +
  scale_shape_manual(values=c(1:16)) +
  theme_minimal()
ggsave("type2_composite.png", width=8, height=5)

#arrange data into table suitable for ggplot, plot, and save (type3) 7
type3.plot.table=c()
type3.plot.table=as.data.frame(
  cbind(
    c(type3$`1G.WT`,type3$MG.WT,type3$`1G.TR`,type3$MG.TR), 
    rep(rownames(type3),4), 
    c(rep("1G.WT",7),rep("MG.WT",7),rep("1G.TR",7),rep("MG.TR",7))
  ), stringsAsFactors=FALSE)
colnames(type3.plot.table)=c("scaled_cts", "gene", "condition")
type3.plot.table$scaled_cts=as.numeric(type3.plot.table$scaled_cts)
type3.plot.table$condition=factor(type3.plot.table$condition, levels = c("1G.WT", "MG.WT", "1G.TR", "MG.TR"))
b = ggplot(type3.plot.table, aes(x=condition, y=scaled_cts, group=gene)) +
  scale_y_continuous(name="Normalized Counts (Scaled)", limits=c(0.4, 1.1)) +
  geom_jitter(position=position_jitter(0.3), aes(shape=gene)) +
  scale_shape_manual(values=c(1:7)) +
  theme_minimal()
ggsave("type3_composite.png", width=8, height=5)

#arrange data into table suitable for ggplot, plot, and save (type4) 5
type4.plot.table=c()
type4.plot.table=as.data.frame(
  cbind(
    c(type4$`1G.WT`,type4$MG.WT,type4$`1G.TR`,type4$MG.TR), 
    rep(rownames(type4),4), 
    c(rep("1G.WT",5),rep("MG.WT",5),rep("1G.TR",5),rep("MG.TR",5))
  ), stringsAsFactors=FALSE)
colnames(type4.plot.table)=c("scaled_cts", "gene", "condition")
type4.plot.table$scaled_cts=as.numeric(type4.plot.table$scaled_cts)
type4.plot.table$condition=factor(type4.plot.table$condition, levels = c("1G.WT", "MG.WT", "1G.TR", "MG.TR"))
b = ggplot(type4.plot.table, aes(x=condition, y=scaled_cts, group=gene)) +
  scale_y_continuous(name="Normalized Counts (Scaled)", limits=c(0.6, 1.1)) +
  geom_jitter(position=position_jitter(0.2), aes(shape=gene)) +
  scale_shape_manual(values=c(1:5)) +
  theme_minimal()
ggsave("type4_composite.png", width=8, height=5)

#arrange data into table suitable for ggplot, plot, and save (type5) 4
type5.plot.table=c()
type5.plot.table=as.data.frame(
  cbind(
    c(type5$`1G.WT`,type5$MG.WT,type5$`1G.TR`,type5$MG.TR), 
    rep(rownames(type5),4), 
    c(rep("1G.WT",4),rep("MG.WT",4),rep("1G.TR",4),rep("MG.TR",4))
  ), stringsAsFactors=FALSE)
colnames(type5.plot.table)=c("scaled_cts", "gene", "condition")
type5.plot.table$scaled_cts=as.numeric(type5.plot.table$scaled_cts)
type5.plot.table$condition=factor(type5.plot.table$condition, levels = c("1G.WT", "MG.WT", "1G.TR", "MG.TR"))
b = ggplot(type5.plot.table, aes(x=condition, y=scaled_cts, group=gene)) +
  scale_y_continuous(name="Normalized Counts (Scaled)", limits=c(0.6, 1.1)) +
  geom_jitter(position=position_jitter(0.25), aes(shape=gene)) +
  scale_shape_manual(values=c(1:4)) +
  theme_minimal()
ggsave("type5_composite.png", width=8, height=5)
