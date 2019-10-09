library(reshape2)
library(Matrix)

####################### HELP function #########################
# read celltype matrix: 
# read table file and set colnames as rownames
read_celltype_mat <- function(cell_type_file){
   message("loading file")

   cell_type_df <- read.table(cell_type_file, sep="\t")
   colnames(cell_type_df) <- rownames(cell_type_df)

   cell_type_mat<-as.matrix(cell_type_df)
   cell_type_mat<-melt(cell_type_mat)

   return(cell_type_mat)
}

# merge matrix:
# rbind
merge_matrix <- rbind

# anno_Cor_file:
# read table-sep annotation information
# cind and annotated as "CellType1", "Species1", "Sub_Cluster1", 
# "CellType2", "Species2", "Sub_Cluster2", "Value"
anno_Cor_file <- function(total_celltype_df, anno_info_fname){
  message("annotating celltype")

  #file IO
  anno_info_df <- read.table(anno_info_fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  row.names(anno_info_df) <- anno_info_df[,1]
  total_celltype_df <- as.data.frame(total_celltype_df)
  total_celltype_df[,1] <- as.character(total_celltype_df[,1])
  total_celltype_df[,2] <- as.character(total_celltype_df[,2])
  
  #check point
  if (CHECK_POINT){
    message("total_celltype_df:")
    print(str(total_celltype_df))
    message("anno_info_df:")
    print(str(anno_info_df))
  }

  #annotation
  anno_celltype1 <- anno_info_df[total_celltype_df[,1], c(1,2,4)]
  anno_celltype2 <- anno_info_df[total_celltype_df[,2], c(1,2,4)]
  anno_celltype_df <- cbind(anno_celltype1, anno_celltype2, total_celltype_df[, 3])

  colnames(anno_celltype_df) <- c("CellType1", "Species1", "Sub_Cluster1", 
                                  "CellType2", "Species2", "Sub_Cluster2", "Value")
  rownames(anno_celltype_df)<-NULL

  anno_celltype_df <- na.omit(anno_celltype_df)
  anno_celltype_df<-anno_celltype_df[anno_celltype_df$Species1 != anno_celltype_df$Species2,]

  #check point
  if (CHECK_POINT){
    message("anno_celltype1:")
    print(str(anno_celltype1))
    message("anno_celltype2:")
    print(str(anno_celltype2))
    message("anno_celltype:")
    print(str(anno_celltype_df))
  }

  #为什么类型转回matrix
  anno_celltype_mat<-as.matrix(anno_celltype_df)
  return(anno_celltype_mat)
}

# celltype_rearrange:
# rearrange celltype using the ordered species_name
# arrange_first: the primer species should be in col2               DEFAULT: FALSE
# the other species should be in col5 and rearranged
celltype_rearrange <- function(anno_cor_mat, species_name_str, arrange_first=FALSE){
	if (arrange_first){
		message("rearrange the first prime celltype")
		n <- which(anno_cor_mat[,2] == species_name_str)
	}else{
		message("rearrange selected celltype")
		n <- which(anno_cor_mat[,5] == species_name_str)
	}
	rearranged_mat <- rbind(anno_cor_mat[-n,], anno_cor_mat[n, c(4,5,6,1,2,3,7)])

	return(rearranged_mat)
}


# delete some celltype
delete_celltype <- function(rearranged_celltype_mat, delete_list){
  rearranged_celltype_mat <- as.data.frame(rearranged_celltype_mat)
  for(item in delete_list){
    rearranged_celltype_mat <- rearranged_celltype_mat[rearranged_celltype_mat$CellType1 != as.character(item),]
    rearranged_celltype_mat <- rearranged_celltype_mat[rearranged_celltype_mat$CellType2 != as.character(item),]
  }
  return(as.matrix(rearranged_celltype_mat))
}

# class_class_similarity_counting:
# counting average_value of Sub_Cluster1 and Sub_Cluster2
# set average_value_celltype_mat > 0.8 as blackground
# counting MAX average_value of sub_cluster1
class_class_similarity_counting <- function(rearranged_celltype_mat){
  #####class_class_similarity
  rearranged_celltype_mat <- as.data.frame(rearranged_celltype_mat, stringsAsFactors=FALSE)
  rearranged_celltype_mat$Value <- as.numeric(rearranged_celltype_mat$Value)

  ##counting average_value
  average_value_celltype_mat <- aggregate(rearranged_celltype_mat[,7],
										by=list(rearranged_celltype_mat[,3], rearranged_celltype_mat[,6]), 
										FUN=mean, na.rm=TRUE)
  colnames(average_value_celltype_mat) <- c("Sub_Cluster1", "Sub_Cluster2", "Value")

  #####MAX_ADD
  ##as the blackground of class simialrity
  ##average_value{sub_cluster1}{sub_cluster2} > 0.8
  blackground_aveVal_celltype <- average_value_celltype_mat[average_value_celltype_mat$Value>0.8, ]

  ##MAX average_value of {sub_cluster1}
  max_aveVal_celltype <- aggregate(average_value_celltype_mat,
                  by=list(average_value_celltype_mat$Sub_Cluster1),
                  FUN=max, na.rm=TRUE, drop=FALSE)
  max_aveVal_celltype <- max_aveVal_celltype[,c(2,3,4)]
  
  ##total celltype anno_info of
  ##average_value{sub_cluster1}{sub_cluster2} > 0.8 OR MAX average_value of {sub_cluster1}
  max_rearranged_celltype <- aggregate(rearranged_celltype_mat,
                    by=list(rearranged_celltype_mat[,3], rearranged_celltype_mat[,6]),
                    FUN=max, na.rm=TRUE, drop=FALSE)
  max_rearranged_celltype <- na.omit(max_rearranged_celltype)
  max_rearranged_celltype <- max_rearranged_celltype[,c(-1,-2)]
  
  #use pasted cluster name as index
  rownames(max_rearranged_celltype) <- paste(max_rearranged_celltype$Sub_Cluster1, max_rearranged_celltype$Sub_Cluster2, sep="")
  blackground_key <- paste(blackground_aveVal_celltype$Sub_Cluster1, blackground_aveVal_celltype$Sub_Cluster2, sep="")
  max_key <- paste(max_aveVal_celltype$Sub_Cluster1, max_aveVal_celltype$Sub_Cluster2, sep="")
  max_black_celltype <- max_rearranged_celltype[union(blackground_key, max_key),]
  rownames(max_black_celltype)<-NULL
  
  if (CHECK_POINT){
    message("average_value_celltype_mat:")
    print(str(average_value_celltype_mat))
    message("blackground_aveVal_celltype")
    print(str(blackground_aveVal_celltype))
    message("max_aveVal_celltype:")
    print(str(max_aveVal_celltype))
    message("max_black_celltype:")
    print(str(max_black_celltype))
  }

  return(list(average_value_celltype_mat, blackground_aveVal_celltype, max_aveVal_celltype, max_black_celltype))
}




########################### Pipeline main function ##########################
# merge_celltype_pipeline:
# fname_vector:   celltype.NV_SRS file, convert to celltype-celtype matrix
# arranged_species_name_vector:   ordered species name, from primer to advance
# CHECK_POINT:    check and write all the temp file,    DEFAULT:  FALSE
merge_celltype_pipeline <- function(fname_vector, arranged_species_name_vector, anno_info_fname, root_fname,
                                    delete_list=NULL, CHECK_POINT=FALSE){
  #if (length(file_name_str) != length(arranged_species_name_vector)){
  #  exit(1)
  #}
  CHECK_POINT <<- CHECK_POINT
  
  total_celltype_mat <- read_celltype_mat(fname_vector[1])
  for (file_name_str in fname_vector[-1]){
    print(file_name_str)
    tmp_celltype_mat <- read_celltype_mat(file_name_str)
    total_celltype_mat <- merge_matrix(total_celltype_mat, tmp_celltype_mat)
  } 
  if (CHECK_POINT) {
      #check point
      write.table(total_celltype_mat, file="Total_dup_species.Cor.txt", 
                  sep="\t", quote=F, row.names=F)
  }

  annotated_total_celltype_mat <- anno_Cor_file(total_celltype_mat, anno_info_fname)
  if (CHECK_POINT) {
    #check point
    write.table(annotated_total_celltype_mat, file="Total_dup_species.Cor.ann_subcluster.txt", 
                sep="\t", quote=F, row.names=F)
  }

  rearranged_anno_total_celltype_mat <- celltype_rearrange(annotated_total_celltype_mat, arranged_species_name_vector[1],
                                                                arrange_first=TRUE)
  if (CHECK_POINT){
    #check point
    write.table(rearranged_anno_total_celltype_mat, file="Total_dup_species.Cor.ann_subcluster1.sort.txt",
                sep="\t", row.names=FALSE, quote=F)
  }
  for (species_name in arranged_species_name_vector[-1]){
    print(species_name)
    rearranged_anno_total_celltype_mat <- celltype_rearrange(rearranged_anno_total_celltype_mat, species_name,
                                                                arrange_first=FALSE)
  }
  if (CHECK_POINT){
      #check point
      write.table(rearranged_anno_total_celltype_mat, file="Total_dup_species.Cor.ann_subcluster.sort.txt",
                  sep="\t", row.names=FALSE, quote=F)
  }

  if (!is.null(delete_list)){
    rearranged_anno_total_celltype_mat <- delete_celltype(rearranged_anno_total_celltype_mat, delete_list)
  }
  ## delete the duplicated rows
  rearranged_anno_total_celltype_mat <- rearranged_anno_total_celltype_mat[duplicated(rearranged_anno_total_celltype_mat),]
  #class_similarity_counting_result <- list(average_value_celltype_mat, blackground_aveVal_celltype, max_aveVal_celltype, max_black_celltype) 
  class_similarity_counting_result <- class_class_similarity_counting(rearranged_anno_total_celltype_mat)
  
  if (!is.null(root_fname)){
    root <- read.table(root_fname, sep='\t')
    colnames(root) <- c("CellType1", "Species1", "Sub_Cluster1", 
                        "CellType2", "Species2", "Sub_Cluster2", "Value")
    class_similarity_counting_result[[4]] <- rbind(class_similarity_counting_result[[4]], root)
  }
  
  if (CHECK_POINT){
      #check point
      write.table(class_similarity_counting_result[[1]], file="Total_dup_species.subclass.Cor.txt",
                  sep="\t", row.names=FALSE, quote=F)
      write.table(class_similarity_counting_result[[2]], file="TT_08_1.out",
                  sep="\t", row.names=FALSE, quote=F)
      write.table(class_similarity_counting_result[[3]], file="TT_08.out",
                  sep="\t", row.names=FALSE, quote=F)
      write.table(class_similarity_counting_result[[4]], file="Total_dup_species.Cor.ann.sort.max_8_subclass.txt",
                  sep="\t", row.names=FALSE, quote=F)
  }else{
    return(class_similarity_counting_result)
  }
}
  


############################ TEST FUNCTION #############################
test_merge_celltype_pipeline <- function(){
  fname_vector <- c("/home/ggj/jiaqiLi/dev/R_dev/JiaQi/H_M/ensemble/AddFat2/celltype.NV_SRS_75.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/M_Z/AddAdipose/M_adult/celltype.NV_SRS_75.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/M_Z/AddAdipose/M_larvea/celltype.NV_SRS_75.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/M_Z/AddAdipose/celltype.NV_SRS_75.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/Z_SS/celltype.NV_SRS_75-total.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/Z_SS/Adult/celltype.NV_SRS_75-total.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/Z_SS/Fetal/celltype.NV_SRS_75-total.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/SS_C/celltype.NV_SRS_75-total.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/C_S/celltype.NV_SRS_75.out",
    "/home/ggj/jiaqiLi/dev/R_dev/JiaQi/S_N/celltype.NV_SRS_75.out")

  arranged_species_name_vector <- c(
    "Nematostella",
    "Schmidtea",
    "Celegans",
    "SeaSquirts",
    "Zebrafish",
    "Mouse",
    "Human"
  )
  
  anno_info_fname = "/home/ggj/jiaqiLi/dev/R_dev/To_Jiaqi/Tree/438celltype-NEW-20190728.annotation"
  
  delete_list = c("TNFalpha cells(H).z", "sst1.1 cells(P).z", "Apelin cells(H).z", "Apelin cells(P).z", 
                  "Testicular cell.m", "Prostate gland cell.m", "Luteal cell.m", 
                  "Mammary gland in lactation.m", "Corneal cells.lz", "Ciona31", "Trophoblast progenitor cell.m", 
                  "Unknown(H).z", "Unknown(B).z")
  
  root_fname = "/home/ggj/jiaqiLi/dev/R_dev/script/root.txt"
  
  # print and save all the tmp file
  merge_celltype_pipeline(fname_vector, arranged_species_name_vector, anno_info_fname, root_fname,
                          delete_list=delete_list, 
                          CHECK_POINT=TRUE)
}