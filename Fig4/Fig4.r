library(dplyr)
library(data.table)
library(igraph)


spType <- "IAC"  # AIS， MIA， IAC

edges <- fread(paste0(spType,"_edge.csv"))
head(edges)
colnames(edges)[1:2] <- c("from","to")


nodes <- fread(paste0(spType,"_nodes.csv"))
head(nodes)
colnames(nodes)[1] <- "id"


# 去除edges中的CT
edges <- edges %>% filter(source.type != "CT") %>% filter(target.type != "CT")
nodes <- nodes %>% filter(id %in% c(edges$from, edges$to))


# 再整理Node的phylum信息,将CT和trans整合进phylum中
for(i in c(1:nrow(nodes))){
  if( is.na(nodes$phylum[i]) ){
    nodes$phylum[i] <- nodes$Type[i]
  }
}

edges <- edges %>% filter(!from %in% c("(k)Bacteria","" )) %>% filter( !to %in% c("(k)Bacteria","" )) 
degrees <- data.frame(table(c(edges$from,edges$to)) )
nodes <- merge(nodes, degrees, by.x = "id", by.y = "Var1") 

nodes <- nodes %>% arrange(desc(Freq)) %>% arrange(Type)
nodes$label <- paste(nodes$Type,
                     c(seq(1,table(nodes$Type)[1],by =1), seq(1, table(nodes$Type)[2], by=1)), 
                     sep = "_") 

source("_myStars.function.r")


library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net

# edges parameters #######
maxEdge = 1
minEdge = 0.1
maxAbsR = 1
minAbsR = 0.3

E(net)$width <- rescale(E(net)$abs.r, to =c(minEdge,maxEdge), from = c(minAbsR,maxAbsR))
E(net)$color <- sapply(E(net)$r,
                       function(x){
                         if(x<=0)alpha("skyblue3",0.3) else alpha("indianred",0.3)
                       })

# vertex parameters #######
V(net)$color <- 
  sapply(V(net)$phylum,
         function(x){
           colors_phylum_DF$colors[which(colors_phylum_DF$phylum == x)] 
         })

max(V(net)$Freq); min(V(net)$Freq)
maxCircle = 12
minCircle = 3
maxDegree = 83
minDegree = 1

V(net)$vex.size <- rescale(V(net)$Freq, to =c(minCircle,maxCircle), from = c(minDegree,maxDegree))

V(net)$vex.label <- 
  sapply(1:length(V(net)),
         function(i) {
           if(V(net)$Freq[i] < 30) {  # V(net)$Freq[i] < 30 for AIS; V(net)$Freq[i] < 10 for MIA; V(net)$Freq[i] < 10 for IAC
             NA
           }else if(V(net)$Type[i] == "genus"){
             V(net)$name[i]
           }else {
             V(net)$label[i]
           }
         })

library("rjson")
# Give the input file name to the function.
result <- fromJSON(file = paste0(spType,"_gephi.Fr.nodes.json"))
result$nodes %>% as.data.frame()

jason.coord <- cbind.data.frame(x = sapply(result$nodes, "[[", "x"),
                                y = sapply(result$nodes, "[[", "y"),
                                id = sapply(result$nodes, "[[", "id"))
jason.coord.layout <- jason.coord[match(names(V(net)),jason.coord$id),1:2] %>% as.matrix

pdf(file = paste0(spType, "_genus.trans_igraph.pdf"), 
    width = 8*sqrt(nrow(nodes)/412) , height = 8*sqrt(nrow(nodes)/412)) 
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)
plot(net, vertex.label=V(net)$vex.label,
     vertex.size=(V(net)$vex.size*sqrt(412/nrow(nodes))), 
     vertex.color = V(net)$color, vertex.frame.color="white",
     #vertex.shape = "circle",
     edge.size = E(net)$width, edge.color = E(net)$color,
     layout=jason.coord.layout, edge.curved=0)

dev.off()
