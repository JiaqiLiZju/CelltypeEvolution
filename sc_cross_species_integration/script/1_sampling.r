id<-read.csv("Fig2_Annotation1_V2.csv",sep=",",head=T)



Inter.id<-id[,c(4,10)]
rownames(Inter.id)<-Inter.id[,1]
colnames(Inter.id)<-c("CellID","Celltype")
#Inter1<-Inter[,as.character(Inter.id$CellID)]

#pheno1<-read.table("pheno.out",head=T,sep="\t")
temp<-table(Inter.id$Celltype)


temp<-as.matrix(temp)
rown<-rownames(temp)
Result1<-matrix(,ncol=2) 
Number_celltype<-length(temp[,1])

###
for(i in 1:Number_celltype){         ################
	location<-which(Inter.id$Celltype==as.character(rown[i]))
	if(temp[i]<100){r<-1:temp[i]}
	else{ 
		#tt<-temp[i]*384/(temp[i]+384)
		r<-sample(1:100,replace=F)
	}
	pl<-as.matrix(Inter.id[location[r],])
	Result1<-rbind(Result1,pl)
	
	}
New.id<-Result1[-1,]
dim(New.id)


data<-readRDS("./Mouse.rds")
Inter<-data
#dd<-Result1[order(Result1[,2]),]
New.data<-data[,as.character(rownames(New.id))]
write.table(New.id,"Human_100_anno.id.txt",sep="\t",row.names=F,quote=F)
#write.table(New.data,"Human_100_data.txt",sep="\t",quote=F)
saveRDS(New.data,"Human_100_data.rds")
