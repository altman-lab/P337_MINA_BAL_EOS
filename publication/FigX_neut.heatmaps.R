#### Packages ####
library(tidyverse)
library(limma)
library(Hmisc)
library(ComplexHeatmap)
select <- dplyr::select

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
attach("results/gene_model_fitting.RData")

gene.key <- dat.BAL.abund.norm.voom$genes %>%
  bind_rows(dat.BE.abund.norm.voom$genes) %>%
  select(geneName, hgnc_symbol) %>%
  distinct()

fdr.cutoff <- 0.1

#### Signif NEUT genes - UP ####
genes_neut_bal_up <- lm_neut_BAL$lm %>%
  filter(variable == "NEUT.pct" & FDR < fdr.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_neut_be_up <- lm_neut_BE$lm %>%
  filter(variable == "NEUT.pct" & FDR < fdr.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif NEUT genes - DOWN ####
genes_neut_bal_dn <- lm_neut_BAL$lm %>%
  filter(variable == "NEUT.pct" & FDR < fdr.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_neut_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < fdr.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif FeNO genes - UP ####
genes_feno_bal_up <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < fdr.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_up <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < fdr.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif FeNO genes - DOWN ####
genes_feno_bal_dn <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < fdr.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < fdr.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Significant intersects ####
#Too few to plot
# genes_both_bal_up <- intersect(genes_neut_bal_up, 
#                                c(genes_feno_bal_up,genes_feno_bal_dn))
# genes_both_be_up <- intersect(genes_neut_be_up, 
#                                c(genes_feno_be_up,genes_feno_be_dn))
# 
# genes_both_bal_dn <- intersect(genes_neut_bal_dn, 
#                                c(genes_feno_bal_up,genes_feno_bal_dn))
# genes_both_be_dn <- intersect(genes_neut_be_dn, 
#                               c(genes_feno_be_up,genes_feno_be_dn))

#### Format heatmap data ####
# genes.OI <- c("IL5RA", "CCR3", "SIGLEC8")
genes.OI <- NULL

dat.BAL.rename <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              select(libID, donorID, visit)) %>% 
  filter(visit == "V5") %>% 
  select(gene, donorID, value) %>% 
  mutate(gene = paste("BAL",gene,sep="_"))

dat.BE.rename <- as.data.frame(dat.BE.abund.norm.voom$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BE.abund.norm.voom$targets %>% 
              select(libID, donorID, visit)) %>% 
  filter(visit == "V5") %>% 
  select(gene, donorID, value) %>% 
  mutate(gene = paste("BE",gene,sep="_"))

#Combine bal and be 
temp <- dat.BAL.rename %>% 
  filter(gene %in% c(genes_neut_bal_up,
                     genes_feno_bal_up,
                     genes_neut_bal_dn,
                     genes_feno_bal_dn)) %>% 
  mutate(NEUT = case_when(gene %in% c(genes_neut_bal_up, genes_neut_bal_dn) ~ "FDR < 0.1",
                          TRUE~"NS"),
         FeNO = case_when(gene %in% c(genes_feno_bal_up,genes_feno_bal_dn) ~ "FDR < 0.1",
                          TRUE~"NS"))

dat.all.rename <- dat.BE.rename %>% 
  filter(gene %in% c(genes_neut_be_up,
                     genes_feno_be_up,
                     genes_neut_be_dn,
                    genes_feno_be_dn)) %>% 
  mutate(NEUT = case_when(gene %in% c(genes_neut_be_up,genes_neut_be_dn) ~ "FDR < 0.1",
                          TRUE~"NS"),
         FeNO = case_when(gene %in% c(genes_feno_be_up,genes_feno_be_dn) ~ "FDR < 0.1",
                          TRUE~"NS")) %>%
  bind_rows(temp) %>% 
  arrange(donorID) %>% 
  pivot_wider(names_from = donorID) 

#### Up genes correlation ####
hm.mat.up <- dat.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(genes_neut_bal_up,genes_neut_be_up,
                     genes_feno_bal_up,genes_feno_be_up))

hm.mat.up.mat <- hm.mat.up %>% 
  select(-FeNO,-NEUT) %>%
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

#Annotations
col.anno <- HeatmapAnnotation(
  `FeNO` = filter(hm.mat.up, grepl("BE_", gene))$FeNO,
  `neutrophil` = filter(hm.mat.up, grepl("BE_", gene))$NEUT,
  col = list(`neutrophil` = c("FDR < 0.1" = "black","NS" = "white"),
             `FeNO` = c("FDR < 0.1" = "grey","NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = FALSE)
row.anno <- rowAnnotation(
  `FeNO` = filter(hm.mat.up, grepl("BAL_", gene))$FeNO,
  `neutrophil` = filter(hm.mat.up, grepl("BAL_", gene))$NEUT,
  col = list(`neutrophil` = c("FDR < 0.1" = "black","NS" = "white"),
             `FeNO` = c("FDR < 0.1" = "grey","NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = TRUE)

#color
library(circlize)
col_fun = colorRamp2(c(-0,0.5,1), c("white", "orange","darkred"))

set.seed(42)
hm.up <- Heatmap(corr.up.R, 
                 #Titles
                 row_title = "BAL genes",
                 column_title = "BE genes",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 3, #column_km = 3, 
                 use_raster = FALSE, 
                 #Gene labels
                 # row_labels = row.lab.up, column_labels = col.lab.up,
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Annotations
                 top_annotation = col.anno,
                 left_annotation = row.anno,
                 #Legend
                 col=col_fun,
                 heatmap_legend_param = list(direction = "horizontal"))

hm.up <- draw(hm.up, heatmap_legend_side = "bottom",
              padding = unit(c(5, 1, 1, 1), "mm"))

#### Down genes correlation ####
hm.mat.dn <- dat.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(genes_neut_bal_dn,genes_neut_be_dn,
                     genes_feno_bal_dn,genes_feno_be_dn))

hm.mat.dn.mat <- hm.mat.dn %>% 
  select(-FeNO,-NEUT) %>%
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

#### Heatmap - UP ####
# rownames only for genes of interest
row.lab.dn <- sapply(strsplit(rownames(corr.dn.R)," "), `[`, 1)
row.lab.dn[!row.lab.dn %in% genes.OI] <- ""
#Same for columns
col.lab.dn <- sapply(strsplit(colnames(corr.dn.R)," "), `[`, 1)
# col.lab.dn[!col.lab.dn %in% genes.OI] <- ""

#Annotations
col.anno <- HeatmapAnnotation(
  `FeNO` = filter(hm.mat.dn, grepl("BE_", gene))$FeNO,
  `neutrophil` = filter(hm.mat.dn, grepl("BE_", gene))$NEUT,
  col = list(`neutrophil` = c("FDR < 0.1" = "black","NS" = "white"),
             `FeNO` = c("FDR < 0.1" = "grey","NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = FALSE)
row.anno <- rowAnnotation(
  `FeNO` = filter(hm.mat.dn, grepl("BAL_", gene))$FeNO,
  `neutrophil` = filter(hm.mat.dn, grepl("BAL_", gene))$NEUT,
  col = list(`neutrophil` = c("FDR < 0.1" = "black","NS" = "white"),
             `FeNO` = c("FDR < 0.1" = "grey","NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = TRUE)

#color
library(circlize)
col_fun2 = colorRamp2(c(-0,0.5,1), c("white", "lightblue","darkblue"))

set.seed(42)
hm.dn <- Heatmap(corr.dn.R, 
                 #Titles
                 row_title = "BAL genes",
                 column_title = "BE genes",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 3, #column_km = 3, 
                 use_raster = FALSE, 
                 #Gene labels
                 # row_labels = row.lab.dn, column_labels = col.lab.dn,
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Annotations
                 top_annotation = col.anno,
                 left_annotation = row.anno,
                 #Legend
                 col=col_fun2,
                 heatmap_legend_param = list(direction = "horizontal"))

hm.dn <- draw(hm.dn, heatmap_legend_side = "bottom",
              padding = unit(c(5, 1, 1, 1), "mm"))

#### Save ####
png(file="publication/FigX_neut.heatmap.up.png", 
    height=(nrow(corr.up.R)/9)+2, 
    width=(ncol(corr.up.R)/3)+1, res=150, units="in")
draw(hm.up, heatmap_legend_side = "bottom", padding = unit(c(5, 1, 1, 1), "mm"))
dev.off()

png(file="publication/FigX_neut.heatmap.dn.png", 
    height=(nrow(corr.dn.R)/9)+2, 
    width=(ncol(corr.dn.R)/3)+1, res=150, units="in")
draw(hm.dn, heatmap_legend_side = "bottom", padding = unit(c(5, 1, 1, 1), "mm"))
dev.off()
