# mystar ------------------------------------
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}


# #试一试下面的颜色好不好看  ---------------------------------
library(scales)
library(ggsci)
show_col(pal_d3("category10")(10))
show_col(pal_lancet("lanonc")(9))
show_col(pal_nejm("default")(8))
show_col(pal_npg("nrc")(10))
show_col(pal_simpsons("springfield")(16))

tmp <-  c(pal_npg("nrc")(10),pal_simpsons("springfield")(16))
Colors.selected <- tmp[!tmp %in% c("#8A9197FF","#197EC0FF", "#F05C3BFF", "#46732EFF","#71D0F5FF", "#370335FF")][1:14]

phyla <- c("Proteobacteria",   "Gemmatimonadetes", "Cyanobacteria", "Firmicutes",  "Actinobacteria",  
           "Bacteroidetes",  "Fusobacteria", "Planctomycetes", "[Thermi]", "Chloroflexi",     
           "Acidobacteria", "Verrucomicrobia",  "OD1","WPS-2","CT", "trans"  )

colors_phylum_DF <-
  cbind.data.frame(phylum = phyla,
                   colors = c(Colors.selected,"white", "gray" ))  # 暗紫色：#370335FF

colors_phylum_DF$phylum_fct <- factor(colors_phylum_DF$phylum, levels = colors_phylum_DF$phylum)
