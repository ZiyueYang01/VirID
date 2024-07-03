
#/usr/bin/Rscript
args <- commandArgs(TRUE)
label_tree <- args[1]
label_csv <- args[2]
label_grey <- c(args[3])
pic_path <- args[4]

library(ggplot2)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(tidyverse)
library(phangorn)


tree <- read.newick(label_tree)
dd <- read.csv(label_csv, header = TRUE)
dd <- dd[order(dd[,1]),]
dd[dd== ""] <- NA

query<-dd %>% filter(is.na(SuperGroup))

node <- c()

for (a in query$Cluster_sseqid){
  q <- which(get.tree(tree)$tip.label == a)
  node <- append(node,q)
}

data <- data.frame(ID=sort(get.tree(tree)$tip.label), Cluster_taxon = dd$Cluster_taxon)

tree<-phangorn::midpoint(tree)

group1 <- split(data$ID,data$Cluster_taxon)


colourCount <-  length(unique(dd$Cluster_taxon))
col<-colorRampPalette(brewer.pal(6, "Set3"))(colourCount)
taxon_list <- unique(dd$Cluster_taxon)
col <- setNames(col,taxon_list[ taxon_list != paste(label_grey) ])
cols <- c(col, setNames("grey",label_grey))

fontSize <- 0.8
ifelse(length(tree[["tip.label"]])>800,fontSize <- 0.3,ifelse(length(tree[["tip.label"]])>500,fontSize <- 0.5,fontSize <- 0.8))

pic <- groupOTU(tree,group1) %>% ggtree(aes(color=group),size=fontSize-0.2)+geom_tiplab(size=fontSize)+
  scale_color_manual(values = alpha(cols,0.8),
                     breaks = c(unique(dd$Cluster_taxon)),
                     labels = c(unique(dd$Cluster_taxon)),
                     name = "Cluster_taxon")
for ( a in node ){
  pic <- pic + geom_highlight(node=as.numeric(a),fill = "red", alpha = 0.8)
}

pdf(paste0(pic_path, ".pdf"))
pic
dev.off()