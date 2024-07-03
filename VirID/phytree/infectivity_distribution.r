#/usr/bin/Rscript
args <- commandArgs(TRUE)
label_csv <- args[1]
pic_path <- args[2]

library(ggplot2)
library(dplyr)
library(tidyverse)

df <-  read.csv(label_csv, header = TRUE)

df_num <- group_by(df,RdRP_super_group)  %>%
  mutate(Virus_of_human_infection =  case_when((Human == 1 ) ~ 1,
                                                    TRUE ~ 0)) %>%
  mutate(Viruses_associated_with_other_vertebrates =  case_when((Human == 0 & Vertebrates == 1) ~ 1,
                                                    TRUE ~ 0)) %>%
  dplyr::summarize(num_human = sum(Virus_of_human_infection),
                   num_vertebrates = sum(Viruses_associated_with_other_vertebrates),
                   num_Plant = sum(Plant)) %>%
  ungroup()  %>%
  pivot_longer(cols = -RdRP_super_group,names_to  = "type",
               values_to  = "num")

df_num [df_num==0]<-NA

df_num_all <- group_by(df,RdRP_super_group)  %>%
  dplyr::summarize(num_all = n()) %>% ungroup()

# 创建堆叠柱状图
p<-ggplot(df_num, mapping = aes(x = RdRP_super_group,y=num)) +
  geom_bar(df_num_all, mapping = aes(x = RdRP_super_group,y=num_all),stat = "identity",fill="#e4e5e3",color="black") +
  geom_text(df_num_all,mapping = aes(x = RdRP_super_group,y=num_all,label=num_all),size=4,hjust=-0.5)+
  geom_bar(aes(fill=type),stat = "identity",color="black") +
  geom_text(aes(fill=type,label=num),size=4,position=position_stack(vjust=0.5))+
  scale_fill_manual(
    name="Infectvity",
    labels = c("Potential human pathogens","Associated with other vertebrates","Associated with plants"),
    breaks = c("num_human","num_vertebrates","num_Plant"),
    values = c(num_human="#c1121f",num_vertebrates="#e9c46a",
                                num_Plant="#90a955"))+
  labs(title = "Distribution of noval or known viral infectivity", x = "Super-Clade", y = "NUM of Contigs") +
  coord_flip() +
  theme_classic()+theme(legend.position = 'bottom')


pdf(paste0(pic_path, ".pdf"))
p
dev.off()