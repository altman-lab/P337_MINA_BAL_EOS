#### Packages ####
library(limma)
library(Hmisc)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
select <- dplyr::select

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
load("publication/signif_genes.RData")

gene.key <- dat.BAL$genes %>%
  bind_rows(dat.BE$genes) %>%
  select(geneName, hgnc_symbol) %>%
  distinct()

#### Format heatmap data ####
# genes.OI <- c("IL5RA", "CCR3", "SIGLEC8")
genes.OI <- NULL

dat.BAL.rename <- as.data.frame(dat.BAL$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BAL$targets %>% 
              select(libID, donorID, visit)) %>% 
  filter(visit == "V5") %>% 
  select(gene, donorID, value) %>% 
  mutate(gene = paste("BAL",gene,sep="_"))

dat.BE.rename <- as.data.frame(dat.BE$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BE$targets %>% 
              select(libID, donorID, visit)) %>% 
  filter(visit == "V5") %>% 
  select(gene, donorID, value) %>% 
  mutate(gene = paste("BE",gene,sep="_"))

#Combine bal and be 
temp <- dat.BAL.rename %>% 
  filter(gene %in% c(genes_all_bal_up,genes_all_bal_dn))

dat.all.rename <- dat.BE.rename %>% 
  filter(gene %in% c(genes_all_be_up,genes_all_be_dn))  %>%
  bind_rows(temp) %>% 
  arrange(donorID) %>% 
  pivot_wider(names_from = donorID) 

#### Up genes correlation ####
hm.mat.up <- dat.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(genes_all_bal_up,genes_all_be_up))

hm.mat.up.mat <- hm.mat.up %>% 
  column_to_rownames("gene")

corr.up <- rcorr(t(hm.mat.up.mat), type = "pearson")

corr.up.R <- as.data.frame(corr.up$r) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #filter BAL vs BE
  filter(grepl("BAL_", rowname)) %>% 
  filter(grepl("BE_", name)) %>% 
  mutate(across(c(rowname, name), ~gsub("BAL_|BE_", "", .))) %>% 
  #Rename columns
  left_join(gene.key, by=c("name"="geneName")) %>% 
  select(-name) %>% 
  pivot_wider(names_from = hgnc_symbol) %>% 
  #rename rows
  left_join(gene.key, by=c("rowname"="geneName")) %>% 
  select(-rowname) %>% 
  #Fix duplicate HGNC
  mutate(hgnc_symbol = case_when(hgnc_symbol=="ADORA3" & SERPINB2 > 0.5 ~ "ADORA3_1",
                                 hgnc_symbol=="ADORA3" & SERPINB2 < 0.5 ~ "ADORA3_2",
                                 TRUE ~ hgnc_symbol)) %>% 
  column_to_rownames("hgnc_symbol") %>% 
  as.matrix()

#### Heatmap - UP ####
# rownames only for genes of interest
row.lab.up <- sapply(strsplit(rownames(corr.up.R)," "), `[`, 1)
row.lab.up[!row.lab.up %in% genes.OI] <- ""
#Same for columns
col.lab.up <- sapply(strsplit(colnames(corr.up.R)," "), `[`, 1)
# col.lab.up[!col.lab.up %in% genes.OI] <- ""

#color
col_fun = colorRamp2(c(-0,0.5,1), c("white", "orange","darkred"))

set.seed(42)
hm.up1 <- Heatmap(corr.up.R, 
                  #Titles
                  row_title = "BAL gene expression",
                  column_title = "BE gene expression",
                  name = "Pearson",
                  row_title_gp = gpar(fontsize = 10),
                  column_title_gp = gpar(fontsize = 10),
                  #Clustering
                  row_km = 4, #column_km = 3,
                  cluster_rows = TRUE,
                  # row_split = 4,
                  #Gene labels
                  row_names_gp = gpar(fontsize = 8), 
                  column_names_gp = gpar(fontsize = 8),
                  #Legend
                  col=col_fun,
                  heatmap_legend_param = list(direction = "horizontal"))
set.seed(42)
hm.up1 <- draw(hm.up1, heatmap_legend_side = "bottom",
               padding = unit(c(5, 1, 1, 1), "mm"))

#Add cluster names
ord <- row_order(hm.up1) %>% plyr::ldply(., data.frame) %>%
  mutate(name = recode(`.id`,
                       # "1"="Cluster 1: Eosinophil (9 genes)",
                       # "2"="Cluster 2: Dendritic (33 genes)",
                       # "3"="Cluster 3: Eosinophil, B, Mast, NK (27 genes)"
                       "1"="Cluster 1: Eosinophil (8 genes)",
                       "2"="Cluster 2: Dendritic and Lymphocytes (20 genes)",
                       "3"="Cluster 3: Basophil (12 genes)",
                       "4"="Cluster 4: Eosinophil and Mast cells (29 genes)"
  )) %>% 
  arrange(`X..i..`) %>% 
  mutate(color=recode(`.id`,
                      "1"="#44AA99",
                      "2"="#DDCC77",
                      "3"="#882255",
                      "4"="#88CCEE"))
ord %>% count(name)

col.vec <- ord$color
names(col.vec) <- ord$name

row.anno <- rowAnnotation(` ` = ord$name,
                          col = list(` ` = col.vec))

set.seed(42)
hm.up <- Heatmap(corr.up.R, 
                 #Titles
                 row_title = "BAL gene expression",
                 column_title = "BE gene expression",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 4, #column_km = 3, 
                 #Gene labels
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Annotations
                 right_annotation = row.anno,
                 #Legend
                 col=col_fun,
                 heatmap_legend_param = list(direction = "horizontal"))

ht_opt$message = FALSE #turn off warning about legend placement
set.seed(42)
hm.up.d <- draw(hm.up, heatmap_legend_side = "bottom",
                padding = unit(c(5, 1, 1, 1), "mm"))

#Save tree for supplemental heatmap
row_dend_up <- row_dend(hm.up.d)
row_ord_up <- row_order(hm.up.d)
#
#### Down genes correlation ####
hm.mat.dn <- dat.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(genes_all_bal_dn,genes_all_be_dn))

hm.mat.dn.mat <- hm.mat.dn %>% 
  column_to_rownames("gene")

corr.dn <- rcorr(t(hm.mat.dn.mat), type = "pearson")

corr.dn.R <- as.data.frame(corr.dn$r) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #filter BAL vs BE
  filter(grepl("BAL_", rowname)) %>% 
  filter(grepl("BE_", name)) %>% 
  mutate(across(c(rowname, name), ~gsub("BAL_|BE_", "", .))) %>% 
  #Rename columns
  left_join(gene.key, by=c("name"="geneName")) %>% 
  select(-name) %>% 
  pivot_wider(names_from = hgnc_symbol) %>% 
  #rename rows
  left_join(gene.key, by=c("rowname"="geneName")) %>% 
  select(-rowname) %>% 
  column_to_rownames("hgnc_symbol") %>% 
  as.matrix()

#### Heatmap - dn ####
# rownames only for genes of interest
row.lab.dn <- sapply(strsplit(rownames(corr.dn.R)," "), `[`, 1)
row.lab.dn[!row.lab.dn %in% genes.OI] <- ""
#Same for columns
col.lab.dn <- sapply(strsplit(colnames(corr.dn.R)," "), `[`, 1)
# col.lab.dn[!col.lab.dn %in% genes.OI] <- ""

#color
col_fun2 = colorRamp2(c(-0,0.5,1), c("white", "skyblue","darkblue"))

set.seed(42)
hm.dn1 <- Heatmap(corr.dn.R,
                 #Titles
                 row_title = "BAL gene expression",
                 column_title = "BE gene Expresssion",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 3, #column_km = 3, 
                 #Gene labels
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Legend
                 col=col_fun2,
                 heatmap_legend_param = list(direction = "horizontal"))
set.seed(42)
hm.dn1 <- draw(hm.dn1, heatmap_legend_side = "bottom",
                padding = unit(c(5, 1, 1, 1), "mm"))

#Add cluster names
ord2 <- row_order(hm.dn1) %>% plyr::ldply(., data.frame) %>%
  mutate(name = recode(`.id`,
                       "1"="Cluster 1: Neutrophil (11 genes)",
                       "2"="Cluster 2: Memory T and Macrophage (15 genes)",
                       "3"="Cluster 3: Macrophage and neutrophil (10 genes)"  )) %>% 
  arrange(`X..i..`) %>% 
  mutate(color=recode(`.id`,
                      "1"="#44AA99",
                      "2"="#DDCC77",
                      "3"="#882255",
                      "4"="#88CCEE"))
ord2 %>% count(name)

col.vec2 <- ord2$color
names(col.vec2) <- ord2$name

row.anno2 <- rowAnnotation(` ` = ord2$name,
                          col = list(` ` = col.vec2))

set.seed(42)
hm.dn <- Heatmap(corr.dn.R, 
                 #Titles
                 row_title = "BAL gene expression",
                 column_title = "BE gene expression",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 3, #column_km = 3, 
                 #Gene labels
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Annotations
                 right_annotation = row.anno2,
                 #Legend
                 col=col_fun2,
                 heatmap_legend_param = list(direction = "horizontal"))

set.seed(42)
hm.dn.d <- draw(hm.dn, heatmap_legend_side = "bottom",
                padding = unit(c(5, 1, 1, 1), "mm"))

#Save tree for supplemental heatmap
row_dend_dn <- row_dend(hm.dn.d)
row_ord_dn <- row_order(hm.dn.d)

#### Save ####

# png(file="publication/Fig2A_eos.heatmap.up.png", 
#     height=10, width=3.5, res=150, units="in")
# hm.up.d
# dev.off()
# 
# png(file="publication/Fig2B_eos.heatmap.dn.png", 
#     height=6, width=4, res=150, units="in")
# hm.dn.d
# dev.off()

pdf(file="publication/Fig2A_eos.heatmap.up.pdf", 
    height=10, width=3.5)
hm.up.d
dev.off()

pdf(file="publication/Fig2B_eos.heatmap.dn.pdf", 
    height=6, width=4)
hm.dn.d
dev.off()

save(row_dend_up, row_dend_dn, 
     row_ord_up, row_ord_dn,
     file="publication/heatmap_dend.RData")

#### Mean corr - UP ####
source("publication/get_hm_cluster.R")
clust_genes_up <- get_hm_cluster(dat = corr.up.R, hm = hm.up.d, 
                              dimension = "row")

corr.up.clust <- as.data.frame(corr.up.R) %>% 
  rownames_to_column("BAL") %>% 
  pivot_longer(-BAL, values_to = "R", names_to = "BE") %>% 
  left_join(clust_genes_up, by=c("BAL"="row")) 

corr.up.clust %>% 
  group_by(cluster) %>% 
  summarise(m = mean(R),
            sd = sd(R))

#Save cluster members for supp table 4
hm_clust_up <- corr.up.clust  %>%
  mutate(BAL = recode(BAL, "ADORA3_1"="ADORA3","ADORA3_2"="ADORA3")) %>% 
  distinct(cluster, BAL)


#### Mean corr - DOWN ####
clust_genes_dn <- get_hm_cluster(dat = corr.dn.R, hm = hm.dn.d, 
                              dimension = "row")

corr.dn.clust <- as.data.frame(corr.dn.R) %>% 
  rownames_to_column("BAL") %>% 
  pivot_longer(-BAL, values_to = "R", names_to = "BE") %>% 
  left_join(clust_genes_dn, by=c("BAL"="row")) 

corr.dn.clust %>% 
  # filter(BE==0."HEY2") %>% 
  group_by(cluster) %>% 
  summarise(m = mean(R),
            sd = sd(R))

#Save cluster members for supp table 4
hm_clust_dn <- corr.dn.clust  %>%
  distinct(cluster, BAL)

#### Save data ####
save(corr.up.R, hm.up.d, corr.dn.R, hm.dn.d,
     hm_clust_up, hm_clust_dn,
     file="publication/heatmap_dat.RData")

