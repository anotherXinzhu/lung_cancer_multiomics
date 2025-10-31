library(data.table)
library(dplyr)
library(reshape2)
library("factoextra")
library("vegan")
library(grid) 

trans_df <- fread("trans_df.fpkm.bm.txt", data.table = F) %>%
  filter(sample != "AAH") %>%
  dcast(sample~Gene, value.var = "bm.value")%>%
  tibble::column_to_rownames("sample")




# grouping information  #######
meta <- fread("ID.match.txt")
meta_df <- meta %>% select(Transcriptome.ID, Group)
meta_df <- meta_df[complete.cases(meta_df),]

all(rownames(trans_df) %in% meta_df$Transcriptome.ID)

min(sapply(meta$Image.ID, nchar),na.rm = T) 
max(sapply(meta$Trans.ID2, nchar),na.rm = T) 
rownames(trans_df) <- 
  sapply(rownames(trans_df),
         function(x){
           if(nchar(x) >= 11){ 
             meta$Sample.ID[grep(x, meta$Image.ID)]
           }else if(nchar(x) <= 7) {
             meta$Sample.ID[which(meta$Trans.ID2 == x  & meta$batch == 2)]
           }else{
             print(x)
           }
         })

if(any(sapply(trans_df,sum) == 0)) dat_positive <- trans_df[,-which(sapply(trans_df,sum) == 0)] else dat_positive <- trans_df
meta_df <-  meta_df[match(rownames(dat_positive), meta_df$Transcriptome.ID ),]

arg.pca <- prcomp(dat_positive, scale. = TRUE) 
summary(arg.pca,loadings = T)

library("factoextra")
get_pca_var(arg.pca)
p <- fviz_eig(arg.pca, addlabels = TRUE)

ind <- get_pca_ind(arg.pca)
ind
p_allvar <- fviz_pca_ind(arg.pca, geom = "point",  pointshape = 21, pointsize = 3,
                         fill = meta_df$Group, 
                         addEllipses = TRUE,
                         legend.title = "sample types",
                         repel = TRUE)
p_allvar


library("vegan")
count = t(dat_positive)
arg.adonis <- adonis2(t(count) ~ Group, data = meta_df, permutations = 999)
arg.adonis$R2[1] ; arg.adonis$`Pr(>F)`[1]

library(grid)
grob <- grobTree(textGrob(paste("adonis: R2 = ",round(arg.adonis$R2[1],3), " p-value = ", arg.adonis$`Pr(>F)`[1], sep = ""), 
                          x=0.4,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=10)))

p1 <- p_allvar + annotation_custom(grob) +
  theme_bw() + theme(panel.grid = element_blank())
p1

