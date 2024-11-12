library(tidyverse)
library(openxlsx)
library(limma)

attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
load("publication/signif_genes.RData")

gene.key <- dat.BAL$genes %>%
  bind_rows(dat.BE$genes) %>%
  select(geneName, hgnc_symbol) %>%
  distinct()

#### Table S1 Demographics ####
meta1 <- dat.BAL$targets %>% 
  distinct(donorID, age_mo, sex, race, FEV1.screening.WLAC) %>% 
  arrange(donorID)

#Blood cells
blood <- read_csv("data_raw/addtl.data/P337_cell.counts.csv") %>% 
  select(ptID, contains("blood")) %>% 
  select(-contains("WBC")) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("tissue","cell","units","SBP-Ag"), sep="_") %>% 
  mutate(`SBP-Ag` = recode(`SBP-Ag`,"V4"="Pre","V5"="Post"),
         `SBP-Ag` = factor(`SBP-Ag`, levels = c("Pre","Post"))) %>% 
  group_by(ptID, `SBP-Ag`) %>% 
  mutate(pct = value/sum(value)*100) %>% 
  ungroup() %>% 
  select(ptID, cell, `SBP-Ag`, pct) %>% 
  mutate(cell=paste("blood",cell, "pct",sep=".")) %>% 
  pivot_wider(names_from = cell, values_from = pct) %>% 
  rename(donorID=ptID)

meta2 <- dat.BAL$targets %>% 
  distinct(donorID, 
           FEV1.pctPP.preAlbuterol_V4, FEV1.pctPP.preAlbuterol_V5,
           FeNO.PreBro_V4, FeNO.PreBro_V5) %>% 
  pivot_longer(-donorID) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider() %>% 
  full_join(
    dat.BAL$targets %>% 
      distinct(donorID, visit, EOS.pct, MONO.pct, 
               PMN.pct, LYM.pct, Epi.pct) %>% 
      rename(BAL.EOS.pct=EOS.pct, BAL.MONO.pct=MONO.pct, 
             BAL.PMN.pct=PMN.pct, BAL.LYM.pct=LYM.pct, BAL.Epi.pct=Epi.pct)
  ) %>% 
  mutate(visit = case_when(visit == "V4" ~ "Pre",
                           visit == "V5" ~ "Post"),
         visit = factor(visit, levels = c("Pre","Post"))) %>% 
  rename(`SBP-Ag`=visit) %>% 
  left_join(blood) %>%
  arrange(donorID, `SBP-Ag`)

write.xlsx(list("Demographics"=meta1,
                "SBP-Ag"=meta2), file = "publication/TableS1_demo.xlsx")

#### Summary stats ####
meta2 %>% 
  pivot_longer(-c(donorID, `SBP-Ag`)) %>% 
  group_by(name, `SBP-Ag`) %>% 
  summarise(m = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm=TRUE)) %>% 
  View()

#### Table S2 Sequencing libraries ####
seq1 <- dat.BAL$targets %>% 
  distinct(group, libID, donorID, visit, batch,
           median_cv_coverage, mapped_reads_w_dups, fastq_total_reads) %>% 
  mutate(visit = case_when(visit == "V4" ~ "Pre",
                           visit == "V5" ~ "Post"),
         visit = factor(visit, levels = c("Pre","Post"))) %>% 
  rename(`SBP-Ag`=visit) %>% 
  arrange(donorID, `SBP-Ag`)

seq2 <- dat.BE$targets %>% 
  distinct(group, libID, donorID, visit, batch,
           median_cv_coverage, mapped_reads_w_dups, fastq_total_reads) %>% 
  mutate(visit = case_when(visit == "V4" ~ "Pre",
                           visit == "V5" ~ "Post"),
         visit = factor(visit, levels = c("Pre","Post"))) %>% 
  rename(`SBP-Ag`=visit) %>% 
  arrange(donorID, `SBP-Ag`)

write.xlsx(list("BAL"=seq1, "BE"=seq2),
           file = "publication/TableS2_sequencing.xlsx")

#### Table S4 LM results ####
load("results/gene_model_fitting.RData")

lm.ls <- list()

for(result in rev(ls(pattern = "lm_"))){
  print(result)
  #Nice name for Excel sheet
  if(grepl("visit", result)){
    name <- gsub("visit","SBP-Ag", result)
  } else {
    name <- result
  }
  name <- strsplit(name, split="_")[[1]]
  name <- paste(name[3], name[2]) 
  
  # group
  if(grepl("BAL", result)){
    g <- "BAL"
  } else {
    g <- "BE"
  }
  
  #results
  temp <- get(result)$lm %>% 
    left_join(dat.BAL$genes, by=c("gene"="geneName")) %>% 
    mutate(group = g) %>% 
    select(group, gene, hgnc_symbol, variable, 
           AveExpr, estimate, pval, FDR) %>% 
    mutate(variable = recode(variable, "visit"="SBP-Ag")) %>% 
    filter(variable != "(Intercept)") %>% 
    arrange(-FDR)
  
  lm.ls[[name]] <- temp
}

load("publication/hm_clusters.RData")

all_clust <- hm_clust_up %>% 
  mutate(association="positive") %>% 
  bind_rows(hm_clust_dn %>% 
              mutate(association="negative"))

lm.ls[["Selected DEGs"]] <- data.frame(gene = c(genes_all_bal_up, genes_all_bal_dn,
                                                genes_all_be_up, genes_all_be_dn),
                                       association = c(rep("positive", length(genes_all_bal_up)),
                                                       rep("negative", length(genes_all_bal_dn)),
                                                       rep("positive", length(genes_all_be_up)),
                                                       rep("negative", length(genes_all_be_dn)))) %>% 
  separate(gene, into=c("group","gene"), sep="_") %>% 
  left_join(gene.key, by=c('gene'='geneName')) %>% 
  left_join(all_clust, by=c('hgnc_symbol'="BAL",
                            "association"="association")) %>% 
  mutate(cluster = case_when(hgnc_symbol=="ADORA3"~"Cluster4",
                             TRUE~cluster)) %>% 
  rename(heatmap_cluster=cluster) %>% 
  select(group, association, gene, hgnc_symbol, heatmap_cluster) %>% 
  arrange(group, association, hgnc_symbol)

#Need to add neut and fev
lm.ls[["delta neut"]] <- NULL
lm.ls[["group fev"]] <- NULL
write.xlsx(lm.ls, file = "publication/TableS4_linear.models.xlsx")

#### Table S3 Counts ####
count.ls <- list()
count.ls[["BAL"]] <- as.data.frame(dat.BAL$E) %>% 
  rownames("gene")
count.ls[["BE"]] <- as.data.frame(dat.BE$E) %>% 
  rownames("gene")

write.xlsx(count.ls, file = "publication/TableS3_voom.log2cpm.xlsx")

#### Table S5 enrichment ####
load("results/DEG_enrichment.RData")

#FDR<0.1 and overlap>1
e1 <- c5_up_signif %>% 
  mutate(category = "C5 GO BP", .before = group) %>% 
  mutate(heatmap = "Postive EOS and FeNO", .before = group) %>% 
  rename(DEG_in_pathway=group_in_pathway,
         DEG_cluster=group) %>% 
  bind_rows(c5_dn_signif %>% 
              mutate(category = "C5 GO BP", .before = group) %>% 
              mutate(heatmap = "Negative EOS and FeNO", .before = group) %>% 
              rename(DEG_in_pathway=group_in_pathway,
                     DEG_cluster=group)) %>% 
  arrange(heatmap, FDR)

e2 <- c7_up_signif %>% 
  mutate(category = "C7 GSE3982", .before = group) %>% 
  mutate(heatmap = "Postive EOS and FeNO", .before = group) %>% 
  rename(DEG_in_pathway=group_in_pathway,
         DEG_cluster=group) %>% 
  bind_rows(c7_up_signif %>% 
              mutate(category = "C7", .before = group) %>% 
              mutate(heatmap = "Negative EOS and FeNO", .before = group) %>% 
              rename(DEG_in_pathway=group_in_pathway,
                     DEG_cluster=group)) %>% 
  arrange(heatmap, FDR)

write.xlsx(list("C5 GO BP"=e1, "C7 GSE3982"=e2),
             file = "publication/TableS5_enrichment_0.1FDR.xlsx")

#### Table S6 mediation ####
load("results/mediation.RData")
head(med_result)

#Linear models
med.ls <- list()

for(out in c("EOS.pct","FeNO_V5")){
  
  main <- med_lm_result %>% 
    filter(outcome==out & term != "(Intercept)" & model == "main") %>% 
    distinct(outcome, model, term, estimate, p.value) %>% 
    mutate(model = recode(model, "main"=paste0(out,"~CDH26"))) %>% 
    rename(variable=term) %>% 
    mutate(cytokine_mediator=NA) %>% 
    select(outcome, cytokine_mediator, everything())
  
  med <- med_lm_result %>% 
    filter(outcome==out & term != "(Intercept)" & model != "main") %>% 
    select(outcome, cyto, model, term, estimate, p.value) %>% 
    rename(cytokine_mediator=cyto, variable=term) %>% 
    mutate(model = recode(model, 
                          "med_out"=paste0(cytokine_mediator,"~CDH26"),
                          "med"=paste0(out,"~CDH26+",cytokine_mediator))) %>% 
    mutate(variable = case_when(variable=="cyto_value"~"cytokine_mediator",
                                TRUE~variable)) %>% 
    arrange(cytokine_mediator)
  
  med.ls[[paste0(gsub("_V5|.pct","",out), "_lm")]] <- bind_rows(main, med) %>% 
    mutate(outcome=gsub("_V5|.pct","",outcome))
  
}

#ACME
temp <- med_result %>% 
  select(outcome, cyto, name, estimate, p) %>% 
  rename(cytokine_mediator=cyto, p.value=p, variable=name) %>% 
  mutate(outcome=gsub("_V5|.pct","",outcome)) %>% 
  arrange(outcome, cytokine_mediator)

med.ls[["EOS_mediation"]] <- temp %>% filter(outcome=="EOS")
med.ls[["FeNO_mediation"]] <- temp %>% filter(outcome=="FeNO")

write.xlsx(med.ls, file = "publication/TableS6_mediation.xlsx")
