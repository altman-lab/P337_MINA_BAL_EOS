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
  filter(gene %in% c(genes_all_be_up0.3,genes_all_be_dn0.3))  %>%
  bind_rows(temp) %>% 
  arrange(donorID) %>% 
  pivot_wider(names_from = donorID) 

#### Up genes correlation ####
hm.mat.up <- dat.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(genes_all_bal_up,genes_all_be_up0.3))

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

#Reorder
load("publication/heatmap_dend.RData")

row_ord_up_key <- c(rep(1, length(row_ord_up$`1`)), 
                    rep(2, length(row_ord_up$`4`)), 
                    rep(3, length(row_ord_up$`3`)), 
                    rep(4, length(row_ord_up$`2`)))
row_ord_up_format <- c(row_ord_up$`1`,row_ord_up$`4`,
                       row_ord_up$`3`,row_ord_up$`2`)
# row_ord_up_format <- unlist(row_ord_up)
corr.up.R <- corr.up.R[row_ord_up_format,]
# corr.up.R <- corr.up.R[c(1:8, 50:69, 38:49, 9:37),]

#### Heatmap - UP ####
#color
col_fun = colorRamp2(c(-0,0.5,1), c("white", "orange","darkred"))

#Add cluster names
# row_order(hm.up1) %>% plyr::ldply(., data.frame) %>%
ord <- data.frame(`.id`=row_ord_up_key) %>% 
  mutate(name = recode(`.id`,
                       "1"="Cluster 1: Eosinophil (8 genes)",
                       "4"="Cluster 2: Dendritic and Lymphocytes (20 genes)",
                       "3"="Cluster 3: Basophil (12 genes)",
                       "2"="Cluster 4: Eosinophil and Mast cells (29 genes)"
  )) %>% 
  mutate(color=recode(`.id`,
                      "1"="#44AA99",
                      "2"="#88CCEE",
                      "3"="#882255",
                      "4"="#DDCC77"))
ord %>% count(name)

col.vec <- ord$color
names(col.vec) <- ord$name

row.anno <- rowAnnotation(` ` = ord$name,
                          col = list(` ` = col.vec))

col_up_lab <- colnames(corr.up.R)
names(col_up_lab) <- colnames(corr.up.R)
col_up_lab[["RUSC1"]] <- "**RUSC1**"
col_up_lab[["CDH26"]] <- "**CDH26**" 
col_up_lab[["SERPINB2"]] <- "**SERPINB2**" 
col_up_lab[["CDC42EP5"]] <- "**CDC42EP5**" 
col_up_lab[["DOK1"]] <- "**DOK1**" 

col_up_lab[["CCL26"]] <- "**--->** CCL26"

set.seed(42)
hm.up <- Heatmap(corr.up.R, 
                 #Titles
                 row_title = "BAL gene expression",
                 column_title = "BE gene expression",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 cluster_rows = FALSE,
                 row_split = row_ord_up_key,
                 #Gene labels
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #bold genes of interest
                 column_labels = gt_render(col_up_lab),
                 #Annotations
                 right_annotation = row.anno,
                 #Legend
                 col=col_fun,
                 heatmap_legend_param = list(direction = "horizontal"))

ht_opt$message = FALSE #turn off warning about legend placement
set.seed(42)
hm.up.d <- draw(hm.up, heatmap_legend_side = "bottom",
                padding = unit(c(5, 1, 1, 1), "mm"))
#
#### Save ####

png(file="publication/FigS2_eos.heatmap.up0.3.png", 
    height=10, width=7, res=150, units="in")
hm.up.d
dev.off()

pdf(file="publication/FigS2_eos.heatmap.up0.3.pdf", 
    height=10, width=7)
hm.up.d
dev.off()
