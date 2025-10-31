library(data.table)

dat <- fread("final_neg.pos.log2dat.csv", data.table = F)

dat <- dat %>%
  tibble::column_to_rownames("ID") %>%
  t() %>% data.frame() 


dat[1:3,1:3]

# 去除outlier 
outliers <- c("Ctrl151","Ctrl51","Ctrl141","Ctrl61")

dat <- dat %>%
  tibble::rownames_to_column("sample") %>%
  filter(!sample %in% outliers) %>%
  tibble::column_to_rownames("sample")


# grouping information  #######
meta <- fread("ID.match.txt")
meta_df <- meta %>% select(Metabolome_ID, Group)
meta_df <- meta_df[complete.cases(meta_df),]

all(rownames(dat) %in% meta_df$Metabolome_ID)

rownames(dat)[!rownames(dat) %in% meta_df$Metabolome_ID]

commonSps <- intersect(rownames(dat), meta_df$Metabolome_ID)
dat <- dat[match(commonSps, rownames(dat)),]
meta_df <- meta_df[match(commonSps, meta_df$Metabolome_ID),]

if(any(sapply(dat,sum) == 0)) dat_positive <- dat[,-which(sapply(dat,sum) == 0)] else dat_positive <- dat


arg.pca <- prcomp(dat_positive, scale. = TRUE) 
summary(arg.pca,loadings = T)

library("factoextra")
get_pca_var(arg.pca)

p <- fviz_eig(arg.pca, addlabels = TRUE)
p

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
