library(tidyverse)
library(limma)
library(Hmisc)
library(ComplexHeatmap)
library(ggpubr)
library(patchwork)
select <- dplyr::select

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
attach("results/gene_model_fitting.RData")
attach("results/signif_gene_lists.RData")

gene.key <- dat.BAL.abund.norm.voom$genes %>% 
  bind_rows(dat.BE.abund.norm.voom$genes) %>% 
  select(geneName, hgnc_symbol) %>% 
  distinct()

#### Stats for text ####
targets.all <- bind_rows(dat.BAL.abund.norm.voom$targets,
                         dat.BE.abund.norm.voom$targets) %>% 
  distinct(donorID, visit, EOS.pct, FeNO.PreBro_V4, FeNO.PreBro_V5) %>% 
  pivot_longer(-c(donorID, visit, EOS.pct), values_to = "FeNO") %>% 
  separate(name, into=c("variable","visit2"), sep="_") %>% 
  filter(visit==visit2) %>% 
  select(donorID, visit, EOS.pct, FeNO) 

targets.all %>% 
  pivot_longer(-c(donorID, visit), names_to = "variable") %>% 
  group_by(visit, variable) %>% 
  summarise(m = mean(value, na.rm=TRUE),
            sd = sd(value, na.rm=TRUE))

targets.all %>% 
  summarise(across(EOS.pct:FeNO,
                   ~t.test(.[visit == "V4"], .[visit == "V5"], paired = TRUE)$p.value))

#### EOS FeNO plot ####
targets.all.format <- targets.all %>% 
  pivot_longer(-c(donorID, visit), names_to = "variable") %>% 
  mutate(fdr = recode(variable, "EOS.pct"="FDR = 4.3E-5",
                           "FeNO"="FDR = 4.1E-5"),
         visit = fct_recode(visit, "Pre"="V4", "Post"="V5")) %>% 
  group_by(donorID, variable) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) 

p1 <- targets.all.format %>% 
  filter(variable == "EOS.pct") %>% 
  ggplot() +
  aes(x=visit, y= value) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  #Beautify
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(linewidth=0.5)) +
  labs(y="Eosinophils (%)", x="SBP-Ag", color="Post - Pre\nchange") +
  facet_grid(~fdr, scales="free") +
  theme(legend.position = "none")
  
p2 <- targets.all.format %>% 
  filter(variable == "FeNO") %>% 
  ggplot() +
  aes(x=visit, y= value) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  #Beautify
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(linewidth=0.5)) +
  labs(y="FeNO", x="SBP-Ag", color="Post - Pre\nchange") +
  facet_grid(~fdr, scales="free")

p_all <- p1+p2
ggsave("publication/abstract/Fig1_EOS.and.FeNO.pdf", p_all, width=5, height=3)

#### Heatmap ####
genes.OI <- c("IL5RA", "CCR3", "SIGLEC8")

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

#Combine 
temp <- dat.BAL.rename %>% 
  filter(gene %in% c(BAL_eos_0.1_up, BAL_eos_0.1_dn)) %>% 
  mutate(EOS = ifelse(gene %in% c(BAL_eos_0.1_up, BAL_eos_0.1_dn), "FDR < 0.1", "NS"),
         FeNO = ifelse(gene %in% c(BAL_feno_0.1_up, BAL_feno_0.1_dn), "FDR < 0.1", "NS"))

dat.visit.all.rename <- dat.BE.rename %>% 
  filter(gene %in% c(BE_eos_0.1_up, BE_eos_0.1_dn)) %>% 
  mutate(EOS = ifelse(gene %in% c(BE_eos_0.1_up, BE_eos_0.1_dn), "FDR < 0.1", "NS"),
         FeNO = ifelse(gene %in% c(BE_feno_0.1_up, BE_feno_0.1_dn), "FDR < 0.1", "NS")) %>% 
  bind_rows(temp) %>% 
  arrange(donorID) %>% 
  pivot_wider(names_from = donorID) 

hm.mat.up <- dat.visit.all.rename %>% 
  select_if(~all(!is.na(.))) %>%
  filter(gene %in% c(BAL_eos_0.1_up,BE_eos_0.1_up))

hm.mat.up.mat <- hm.mat.up %>% 
  select(-EOS, -FeNO) %>% 
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
  mutate(hgnc_symbol = case_when(hgnc_symbol=="ADORA3" & AK2 > 0.64 ~ "ADORA3_1",
                                 hgnc_symbol=="ADORA3" & AK2 < 0.64 ~ "ADORA3_2",
                                 TRUE ~ hgnc_symbol)) %>% 
  column_to_rownames("hgnc_symbol") %>% 
  as.matrix()

col.anno <- HeatmapAnnotation(
  `BE EOS FDR < 0.1` = filter(hm.mat.up, grepl("BE_", gene))$EOS,
  `BE FeNO FDR < 0.1` = filter(hm.mat.up, grepl("BE_", gene))$FeNO,
  col = list(`BE EOS FDR < 0.1` = c("FDR < 0.1" = "#E1BE6A", "NS" = "white"),
             `BE FeNO FDR < 0.1` = c("FDR < 0.1" = "#40B0A6", "NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = c(FALSE,FALSE))
row.anno <- rowAnnotation(
  `BAL EOS FDR < 0.1` = filter(hm.mat.up, grepl("BAL_", gene))$EOS,
  `BAL FeNO FDR < 0.1` = filter(hm.mat.up, grepl("BAL_", gene))$FeNO,
  col = list(`BAL EOS FDR < 0.1` = c("FDR < 0.1" = "#E1BE6A", "NS" = "white"),
             `BAL FeNO FDR < 0.1` = c("FDR < 0.1" = "#40B0A6", "NS" = "white")),
  annotation_name_gp= gpar(fontsize = 8),
  simple_anno_size = unit(0.2, "cm"), 
  show_legend = c(FALSE,FALSE))

# rownames only for genes of interest
row.lab <- sapply(strsplit(rownames(corr.up.R)," "), `[`, 1)
row.lab[!row.lab %in% genes.OI] <- ""
#Same for columns
col.lab <- sapply(strsplit(colnames(corr.up.R)," "), `[`, 1)
# col.lab[!col.lab %in% genes.OI] <- ""

#color
library(circlize)
col_fun = colorRamp2(c(0,0.5,1), c("white", "orange","darkred"))

set.seed(42)
hm.up <- Heatmap(corr.up.R, 
                 #Titles
                 row_title = "BAL genes associated with eosinophils",
                 column_title = "BE genes associated with eosinophils",
                 name = "Pearson",
                 row_title_gp = gpar(fontsize = 10),
                 column_title_gp = gpar(fontsize = 10),
                 #Clustering
                 row_km = 3, column_km = 3, 
                 use_raster = FALSE, 
                 #Gene labels
                 row_labels = row.lab, column_labels = col.lab,
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8),
                 #Annotations
                 top_annotation = col.anno, left_annotation = row.anno,
                 
                 #Legend
                 col=col_fun,
                 heatmap_legend_param = list(direction = "horizontal"))

pdf(file="publication/abstract/Fig2_EOSup.heatmap.pdf", height=8, width=6)
draw(hm.up, heatmap_legend_side = "bottom", padding = unit(c(5, 1, 1, 1), "mm"))
dev.off()

#### Gene plots: expression pre/post, EOS, FeNO ####
genes.OI2 <- c("SERPINB2", "CCR3")

dat.feno <- targets.all

#Combine all expression data
dat.bal <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  #get gene symbols
  rownames_to_column("geneName") %>% 
  left_join(select(dat.BAL.abund.norm.voom$genes, geneName, hgnc_symbol)) %>% 
  #add and format metadata
  pivot_longer(-c(geneName, hgnc_symbol), names_to = "libID") %>% 
  left_join(select(dat.BAL.abund.norm.voom$targets, libID, donorID, visit)) %>% 
  left_join(targets.all) %>% 
  mutate(visit = fct_recode(visit, "Pre"="V4", "Post"="V5"),
         cell = "BAL") %>%
  select(geneName, hgnc_symbol, donorID, libID, cell, visit, EOS.pct, FeNO, value) %>% 
  #genes of interest
  filter(hgnc_symbol %in% genes.OI2)

dat.be <- as.data.frame(dat.BE.abund.norm.voom$E) %>% 
  #get gene symbols
  rownames_to_column("geneName") %>% 
  left_join(select(dat.BE.abund.norm.voom$genes, geneName, hgnc_symbol)) %>% 
  #add and format metadata
  pivot_longer(-c(geneName, hgnc_symbol), names_to = "libID") %>% 
  left_join(select(dat.BE.abund.norm.voom$targets, libID, donorID, visit)) %>% 
  left_join(targets.all) %>% 
  mutate(visit = fct_recode(visit, "Pre"="V4", "Post"="V5"),
         cell = "BE") %>%
  select(geneName, hgnc_symbol, donorID, libID, cell, visit, EOS.pct, FeNO, value) %>% 
  #genes of interest
  filter(hgnc_symbol %in% genes.OI2)

#Format lme results
fdr.bal <- bind_rows(lm_visit_BAL$lm, lm_eos_BAL$lm, lm_feno_BAL$lm) %>% 
  filter(variable %in% c("visit", "EOS.pct", "FeNO") & gene %in% dat.bal$geneName) %>% 
  #get gene symbols
  left_join(select(dat.BAL.abund.norm.voom$genes, geneName, hgnc_symbol), 
            by=c("gene"="geneName")) %>%
  select(hgnc_symbol, variable, estimate, FDR) %>%
  pivot_longer(estimate:FDR) %>% 
  mutate(name = paste(variable, name, sep="_")) %>% 
  select(-variable) %>% 
  pivot_wider() %>% 
  #Create start/stop groups for FDR bars
  mutate(cell = "BAL") 

fdr.be <- bind_rows(lm_visit_BE$lm, lm_eos_BE$lm, lm_feno_BE$lm) %>% 
  filter(variable %in% c("visit", "EOS.pct", "FeNO") & gene %in% dat.bal$geneName) %>% 
  #get gene symbols
  left_join(select(dat.BE.abund.norm.voom$genes, geneName, hgnc_symbol), 
            by=c("gene"="geneName")) %>%
  select(hgnc_symbol, variable, estimate, FDR) %>%
  pivot_longer(estimate:FDR) %>% 
  mutate(name = paste(variable, name, sep="_")) %>% 
  select(-variable) %>% 
  pivot_wider() %>% 
  #Create start/stop groups for FDR bars
  mutate(cell = "BE") 

#Add FDR to expression data
dat.all2 <- bind_rows(dat.bal, dat.be) %>% 
  left_join(bind_rows(fdr.bal, fdr.be)) %>% 
  mutate(facet.visit = paste0(cell, "\nFDR = ", 
                              signif(visit_FDR, digits=2)),
         facet.eos = cell,
         facet.feno = cell)

plot.ls <- list()
for(g in genes.OI2){
  plot.temp <- list()
  #Visit
  plot.temp[[paste0(g,"_1")]] <- dat.all2 %>% 
    filter(hgnc_symbol==g) %>%
    ggplot(aes(x=visit, y=value)) +
    geom_path(aes(group=donorID)) +
    geom_point() +
    #Beautify
    theme_classic(base_size = 12) +
    theme(strip.background = element_rect(linewidth=0.5)) +
    labs(y="Normalized log2 expression", x="SBP-Ag", color="",
         title=paste0(LETTERS[which(genes.OI2 == g)],")  ",g)) +
    facet_grid(~facet.visit, scales="free")
  #EOS
  plot.temp[[paste0(g,"_2")]] <- dat.all2 %>% 
    filter(visit == "Post") %>% 
    filter(hgnc_symbol==g) %>%
    ggplot(aes(x=EOS.pct, y=value)) +
    geom_point(color="grey80") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = 'y ~ x')+
    #Beautify
    theme_classic(base_size = 12) +
    theme(strip.background = element_rect(linewidth=0.5)) +
    labs(y="Normalized log2 expression", x="Eosinophil (%)", color="") +
    facet_grid(~facet.eos, scales="free") +
    stat_cor(method = "pearson", aes(label = ..r.label..),
             label.x.npc = "right", label.y.npc = "bottom", hjust=1)
  
  #FeNO
  plot.temp[[paste0(g,"_3")]] <- dat.all2 %>% 
    filter(visit == "Post") %>% 
    filter(hgnc_symbol==g) %>%
    ggplot(aes(x=FeNO, y=value)) +
    geom_point(color="grey80") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = 'y ~ x')+
    #Beautify
    theme_classic(base_size = 12) +
    theme(strip.background = element_rect(linewidth=0.5)) +
    labs(y="Normalized log2 expression", x="FeNO", color="") +
    facet_grid(~facet.feno, scales="free") +
    stat_cor(method = "pearson", aes(label = ..r.label..),
             label.x.npc = "right", label.y.npc = "bottom", hjust=1)
  
  plot.ls[[g]] <- wrap_plots(plot.temp) +plot_layout(ncol = 1)
}

wrap_plots(plot.ls)
ggsave("publication/abstract/Fig3_genes.pdf", wrap_plots(plot.ls), 
       height=7.5, width=6)
