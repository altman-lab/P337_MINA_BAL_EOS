library(tidyverse)

load("publication/signif_genes.RData")

#### Count concordant DEGs ####
gene.ls <- list()
gene.ls[["BAL_all_dn"]] <- genes_all_bal_dn
gene.ls[["BAL_all_up"]] <- genes_all_bal_up
gene.ls[["BE_all_dn"]] <- genes_all_be_dn
gene.ls[["BE_all_up"]] <- genes_all_be_up

gene.df <- plyr::ldply(gene.ls, data.frame) %>% 
  rename(group=`.id`, ID=`X..i..`) %>% 
  separate(group, into = c("cell","variable","direction"), sep="_") %>% 
  separate(ID, into = c("cell2","gene"), sep="_", remove = FALSE) %>% 
  pivot_wider(names_from = variable, values_from = direction)

gene.df %>%
  drop_na(all) %>%
  count(cell, all)
