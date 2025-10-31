
# PCA ===========================================
library(data.table)
library(dplyr)
library(reshape2)
library("factoextra")
library("vegan")
library(grid)

otu_df <- fread("16S_L6.bm.txt") %>%
  reshape2::dcast(sample~OTU, value.var = "bm.value") %>%
  tibble::column_to_rownames("sample")

if(any(sapply(otu_df,sum) == 0)) dat_positive <- otu_df[,-which(sapply(otu_df,sum) == 0)] else dat_positive <- otu_df

Meta <- fread("ID.match.txt")

rownames(dat_positive)[!(rownames(dat_positive) %in% Meta$`16SrDNA.ID`) ] <- 
  sapply(rownames(dat_positive)[!(rownames(dat_positive) %in% Meta$`16SrDNA.ID`) ],
         function(x) Meta$`16SrDNA.ID`[grep(x,Meta$Image.ID, fixed = T)])

all(rownames(dat_positive) %in% Meta$`16SrDNA.ID`)

meta_df <- Meta %>% select(`16SrDNA.ID`, Group)
meta_df <- meta_df[complete.cases(meta_df),]

outliers <- c("1A1","137A1","61A1", "67A1")
dat_positive <- dat_positive %>% 
  tibble::rownames_to_column("sample") %>%
  filter(!sample %in% outliers) %>%
  tibble::column_to_rownames("sample")

meta_df <-  meta_df[match(rownames(dat_positive), meta_df$`16SrDNA.ID` ),]

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


# alpha diversity: ===========================================
library(ggpubr) 
library(rlist)

dat <- fread("alphadiversity.txt",select = c(1:2))  
tmp <- t(combn(unique(dat$Group), 2)) 
my_comparisons <- list()
for(i in c(1:nrow(tmp)) ){
  my_comparisons <- list.append(my_comparisons, tmp[i,])
}

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1,1), symbols = c("****", "***", "**", "*", "*","ns"))

dat$Group <- factor(dat$Group, levels = c("BENIGN","AIS","MIA","IAC"))

ggboxplot(dat, x = "Group", y = "shannon",
          color = "Group",width = 0.6, 
          outlier.shape = NA, 
          add="jitter") + 
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     vjust = 0.5,
                     hide.ns=T) +
  theme(legend.position="none",
        axis.text.x=element_text(angle = 50,hjust = 1)) 

# L6 composition in groups : ===========================================
remove(list = ls())
otu_df <- fread("16S_L6.bm.txt") 
Meta <- fread("ID.match.txt")
otu_df$sample[!(otu_df$sample %in% Meta$`16SrDNA.ID`) ] <- 
  sapply(otu_df$sample[!(otu_df$sample %in% Meta$`16SrDNA.ID`) ],
         function(x) Meta$`16SrDNA.ID`[grep(x,Meta$Image.ID, fixed = T)])
all(otu_df$sample %in% Meta$`16SrDNA.ID`)

meta_df <- Meta %>% select(`16SrDNA.ID`, Group)
meta_df <- meta_df[complete.cases(meta_df),]

# define other otu
otuRank_df <- otu_df %>% 
  group_by(OTU) %>%
  summarise(meanRelAbund = mean(relA.aftrBm)) %>%
  arrange(desc(meanRelAbund))

otherOTU <- otuRank_df$OTU[15:nrow(otuRank_df)]  

otu_df$OTU[otu_df$OTU %in% otherOTU] <- "other"

otu_df <- otu_df %>% 
  group_by(sample, OTU) %>%
  summarise(relAbund = sum(relA.aftrBm))

otu_df$Group <-
  sapply(otu_df$sample,
         function(x) meta_df$Group[meta_df$`16SrDNA.ID` == x])

otu.group_df <- otu_df %>% 
  group_by(OTU, Group) %>%
  summarise(mRelAbund = mean(relAbund))

# plotting
library(ggplot2)
unique(otu.group_df$Group)
otu.group_df <- otu.group_df %>% filter(Group != "AAH")

otu.group_df$Group <- 
  factor(otu.group_df$Group, levels = c("BENIGN","AIS","MIA","IAC"))

genus.color_df <- cbind.data.frame(
  genus = c("(g)Acinetobacter","(g)Enhydrobacter","(g)Bacillus","(g)Pseudomonas","(g)Staphylococcus",  
            "(g)Deinococcus","(g)Sphingomonas","(g)Aerococcus","(g)Methylobacterium","(g)Stenotrophomonas",
            "(g)Lactobacillus","(g)Arcobacter","(g)Paracoccus","(g)Micrococcus","other"),
  color = c("#FFFACD","#FF1493","#FFC0CB","#FFA500","#FFFF00",
            "#C0FF3E", "#00CD00", "#008B45", "#87CEFA", "#1E90FF",
            "#000080", "#B452CD","#7A378B","#C1CDC1","gray"),
  stringsAsFactors=F
)


Colors <- sapply(unique(otu.group_df$OTU),
                 function(x)genus.color_df$color[which(genus.color_df$genus == x)])

ggplot(otu.group_df) +
  geom_col(aes(x=Group,y=mRelAbund, fill=OTU),width = 0.6)+
  scale_fill_manual(values = Colors) + 
  theme_bw() + theme(panel.grid = element_blank())

# calculate significance 
GrpPairs <- list(c("BENIGN","AIS"),c("AIS","MIA"), c("MIA","IAC"))

for(gp in GrpPairs){
  # print(gp)
  dat.tmp <- otu_df %>% filter(Group %in% gp)
  
  res <- dat.tmp %>%
    group_by(OTU) %>%
    do(w=wilcox.test(relAbund~Group, data = .)) %>%
    summarise(OTU,wilcoxP = w$p.value) %>%
    arrange(wilcoxP)
  
  writeLines(paste(gp, collapse = " vs. "))
  
  writeLines( paste(res[res$wilcoxP < 0.1,]$OTU, 
                    res[res$wilcoxP < 0.1,]$wilcoxP,
                    sep = "; " ) )
}

# Taxon composition in individual samples: ==============================================
library(reshape2)
library(scales)
dat<-fread("genuslevelAIS_MIA_IAC_BENIGN.txt",header = T,check.names=F,sep = "\t")
ID.match_df<-fread("ID.match.txt",header = T,check.names=F,sep = "\t")

dat.l <- dat %>% reshape2::melt(id.vars="genus", variable.name='sample')
unique(dat.l$genus) 

genus_orderDf <-
  dat.l %>%
  group_by(genus) %>%
  summarize(avgRelAbund = mean(value)) %>%
  arrange(desc(avgRelAbund)) %>%
  mutate(avgRelAbund_p = percent(avgRelAbund, 1))

sum(genus_orderDf$avgRelAbund[15:nrow(genus_orderDf)]) 
othergenus <- genus_orderDf$genus[15:nrow(genus_orderDf)] 
genus_inRank <- c(genus_orderDf$genus[1:14], "other") 

dat.l$genus[dat.l$genus %in% othergenus] <- "other"
dat.l <- dat.l %>% group_by(sample, genus) %>% summarise(relAbund=sum(value))

dat.w <- dat.l %>% reshape2::dcast(sample~genus, value.var = "relAbund")
row.names(dat.w) <- dat.w$sample
dat.w <- dat.w[,-1]

library(ggdendro)

if(T){
  df <- t(dat.w) #自己的数据套到模板里
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x)))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x))))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] #wide dataframe 的行和列重新排列???
  xx_names <- attr(xx, "dimnames") #Clust的顺???
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) #clust图的data
  ddata_y <- dendro_data(dd.col) #clust图的data
} #模板用if(T)包起???

xx_names[[1]]
xx_names[[2]]

dat.l$sample <- factor(dat.l$sample, levels = xx_names[[2]])

library(ggplot2)
library(ggsci)

genus.color_df <- cbind.data.frame(
  genus = c("(g)Acinetobacter","(g)Enhydrobacter","(g)Bacillus","(g)Pseudomonas","(g)Staphylococcus",  
            "(g)Deinococcus","(g)Sphingomonas","(g)Aerococcus","(g)Methylobacterium","(g)Stenotrophomonas",
            "(g)Lactobacillus","(g)Arcobacter","(g)Paracoccus","(g)Micrococcus","other"),
  color = c("#FFFACD","#FF1493","#FFC0CB","#FFA500","#FFFF00",
            "#C0FF3E", "#00CD00", "#008B45", "#87CEFA", "#1E90FF",
            "#000080", "#B452CD","#7A378B","#C1CDC1","gray"),
  stringsAsFactors=F
)

dat.l$genus <- factor(dat.l$genus, levels = genus_inRank)
Colors <- sapply(genus_inRank, 
                 function(x) genus.color_df$color[which(genus.color_df$genus == x)])

ggplot(dat.l, aes(x = sample, y = relAbund, fill = genus)) + 
  geom_col() +  
  scale_fill_manual(values = Colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  #facet_grid(cols=vars(sampleTypes),scales = "free_x", space = "free_x") + 
  xlab("") + ylab("percentage of drug types") + guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position="bottom",
        legend.key.size = unit(0.4, "inches"),
       # legend.title = theme_text(size=50)), 
       panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_none + theme(axis.title.x=element_blank())
