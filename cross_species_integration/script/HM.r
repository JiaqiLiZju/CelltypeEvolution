data1<-readRDS("./Human/Human.rds") ###human dge
raw.10e6<-t(t(data1)/colSums(data1))*1000000
human<-as.data.frame(raw.10e6)

data1<-readRDS("./Mouse/Fig2_Mouse_pse20.rds")  ###mouse dge 
#raw1[raw1<0]=0
raw.10e6<-t(t(data1)/colSums(data1))*1000000
mouse<-as.data.frame(raw.10e6)
#data2<-as.data.frame(data2)

rm(data1, raw.10e6)

orth<-read.table("./Human_Mouse_one-one.ensemble.orth",sep="\t")
orth<-as.data.frame(orth)
mouse.orth<-mouse[orth[,4],]
human.orth<-human[as.character(orth[,2]),]


data<-cbind(mouse.orth,human.orth)
rownames(data)<-orth[,3]
data[is.na(data)]<-0

P1<-read.table("./Phenotype.out",sep="\t",head=T)
data1<-data[,as.character(P1$Sample_ID)]
source("./2017-08-28-runMN-US.R")

celltypes1 <-unique(as.character(P1$Celltype))


var.genes1=get_variable_genes(data1,P1)
length(var.genes1)
write.table(var.genes1,"var.genes_75.out",sep="\t",quote=F)




#####remove------------------------------------------
library(gplots)
library(RColorBrewer)

celltype.NV=run_MetaNeighbor_US(var.genes1,data1,celltypes1,P1)
write.table(celltype.NV,file="celltype.NV_SRS_75.out",sep="\t",quote=F)###---------

cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
pdf("celltype.NV_SRS_75.pdf")  #########--------------------------
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.3,cexCol=0.3)
dev.off()
top_hits=get_top_hits(celltype.NV,P1,threshold=0.9,filename="top_hits_SRS_75.out") 
top_hits=get_top_hits(celltype.NV,P1,threshold=0.8,filename="top_hits_SRS_0.8_75.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.7,filename="top_hits_SRS_0.7_75.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.6,filename="top_hits_SRS_0.6_75.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0,filename="top_hits_SRS_0_75.out")













