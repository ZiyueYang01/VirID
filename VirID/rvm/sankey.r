args <- commandArgs(TRUE)
data <- args[1]
pic_path <- args[2]

library(jsonlite)
library(networkD3)
library(dplyr)
library(ggplot2)

# data <- read.csv(args[1])

nodes <- data.frame(name = unique(c(data$`Node1`, data$`Node2`)))
nodes_group <- unique(merge(nodes, data[ ,colnames(data) %in% c("Node1","group")], by.x = "name", by.y = "Node1",all = TRUE))
nodes_group[is.na(nodes_group)] <-  6
nodes_group<- nodes_group[match(nodes$name, nodes_group$name), ]     


data <- data %>%
  mutate(source = match(`Node1`, nodes$name) - 1,
         target = match(`Node2`, nodes$name) - 1)


links <- data %>%
  group_by(source, target) %>%
  summarise(value = sum(Value)) %>%
  ungroup()

sankey <- sankeyNetwork(Links = links, Nodes = nodes_group, Source = "source", 
                        Target = "target", Value = "value",numberFormat=".1f",
                        fontSize = 12, nodeWidth = 30,NodeGroup="group",showNodeValues=FALSE,fontFamily = "Arial",
                        orderByPath =TRUE,nodeShadow=TRUE,curvature = 0.1,zoom = TRUE,
                        colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"))


htmlwidgets::saveWidget(sankey, pic_path)
