library(networkD3)
library(igraph)
library(tidyr)
library(dplyr)
library(RColorBrewer)

draw_diagonal <- function(simi_cor_data, anno_info_fname, color, save=TRUE, save_name='diag_'){
  data = read.table(simi_cor_data,
                    header=TRUE, sep='\t')
  # check
  #dim(data)
  #colnames(data)
  
  Relationships<- data.frame(Parent=data$CellType2,
                             Child=data$CellType1)
  root<-setdiff(Relationships$Parent, Relationships$Child)
  
  # check
  #root
  #dim(Relationships)
  
  net <- graph_from_data_frame(d=Relationships, directed=T)
  
  # anno
  Clu<-read.table(anno_info_fname,sep="\t",head=T)
  rownames(Clu)<-Clu$Celltype
  Clu$Cluster<-gsub("[.].","",Clu$Cluster)
  Clu$Cluster<-gsub("[..].","",Clu$Cluster)
  
  #traverse next layer and then recurve
  as.list.igraph <- function(thisNode) {
    nm <- vertex_attr(net, "name", thisNode)
    childNodes <- V(net)[which(shortest.paths(net, thisNode, mode="out") == 1)]
    if (length(childNodes)==0) return(list(name=nm))
    list(name=nm, children=unname(lapply(childNodes, as.list.igraph)))
  }
  
  #color
  color.list<-read.table(color,sep="\t",head=T)
  rownames(color.list)<-color.list$name
  Cluster_col<-color.list[Clu$Cluster,]
  color.df<-data.frame(
    name=Clu$Celltype,
    color=Cluster_col$color
  )
  
  w <- paste('{', paste(color.df %>% 
                          mutate(name = paste0('"', name), color = paste0(color, '"')) %>%
                          unite('x', c(name, color), sep = '" : "' ) %>%
                          .$x, collapse= ', '), '}', collapse = '')
  node.col.func <- JS(paste0('function(d, i) { return ', w, '[d.data.name]; }'))
  
  #radialNetwork
  radialNetwork(as.list.igraph(V(net)[root[1]]),
                fontSize =0,
                opacity = 1,
                height = 1000,
                width=1000,
                linkColour = "darkgray",
                textColour = "#cccccc",
                #textColour =node.col.func,
                #nodeStroke = node.col.func,
                nodeColour  = node.col.func)
}


##### TEST Function ####
test_draw_diagonal <- function(){
  draw_diagonal('/home/ggj/jiaqiLi/dev/R_dev/JiaQi/Tree_new_work/Total_dup_species.Cor.ann.sort.max_8_subclass.txt',
                "438celltype-NEW-20190728.annotation",
                "39Cluster_color.major.list")
  }