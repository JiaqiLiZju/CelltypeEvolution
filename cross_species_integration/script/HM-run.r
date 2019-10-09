#!/usr/bin/R

orth<-read.table("./Human_Mouse_one-one.ensemble.orth",sep="\t", stringsAsFactors=F)
rownames(orth) <- orth[,1]
var.genes <- read.table("./var_genes_75_filter.out", sep="\t", stringsAsFactors=F)
var.genes <- as.data.frame(var.genes)
#rownames(var.genes)<-as.data.frame(var.genes)[,1]
orth.var <- orth[var.genes[,1],]
write.table(orth.var, file="var_genes/t3_vargenes.csv", row.names=F, col.names=F, quote=F)

human<-readRDS("./Human/Human_100_data.rds")
###human dge
#raw.10e6<-t(t(data1)/colSums(data1))*1000000
#human<-as.data.frame(raw.10e6)

mouse<-readRDS("./Mouse/Mouse_100_data.rds")
###mouse dge 
#raw1[raw1<0]=0
#raw.10e6<-t(t(data1)/colSums(data1))*1000000
#mouse<-as.data.frame(raw.10e6)
#data2<-as.data.frame(data2)
#rm(data1, raw.10e6)

mouse.orth<-mouse[orth.var[,4],]
mouse.orth[is.na(mouse.orth)]<-0

human.orth<-human[orth.var[,2],]
human.orth[is.na(human.orth)]<-0

rownames(mouse.orth)<-orth.var[,2]
rownames(human.orth)<-orth.var[,2]

write.table(human.orth, file="./Human/Human_vargenes_t3.csv", sep=",", row.names=T, col.names=T, quote=F)
saveRDS(object=human.orth, file="./Human/Human_vargenes_t4.csv.rds")

write.table(mouse.orth, file="./Mouse/Mouse_cell100_vargenes.csv", sep=",", row.names=T, col.names=T, quote=F)
saveRDS(object=mouse.orth, file="./Mouse/Mouse_cell100_vargenes.csv.rds")
