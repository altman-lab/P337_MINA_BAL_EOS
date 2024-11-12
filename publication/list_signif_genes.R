#### Packages ####
library(limma)
library(tidyverse)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
attach("results/gene_model_fitting.RData")

visit.cutoff <- 0.3
bal.cutoff <- 0.1
be.cutoff <- 0.1
be.cutoff2 <- 0.3

#### Signif SBP-Ag genes - UP ####
genes_sbp_bal_up <- lm_visit_BAL$lm %>%
  filter(variable == "visit" & FDR < visit.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_sbp_be_up <- lm_visit_BE$lm %>%
  filter(variable == "visit" & FDR < visit.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif SBP-Ag genes - DOWN ####
genes_sbp_bal_dn <- lm_visit_BAL$lm %>%
  filter(variable == "visit" & FDR < visit.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_sbp_be_dn <- lm_visit_BE$lm %>%
  filter(variable == "visit" & FDR < visit.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif EOS genes - UP ####
genes_eos_bal_up <- lm_eos_BAL$lm %>%
  filter(variable == "EOS.pct" & FDR < bal.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_eos_be_up <- lm_eos_BE$lm %>%
  filter(variable == "EOS.pct" & FDR < be.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

genes_eos_be_up0.3 <- lm_eos_BE$lm %>%
  filter(variable == "EOS.pct" & FDR < be.cutoff2 & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif EOS genes - DOWN ####
genes_eos_bal_dn <- lm_eos_BAL$lm %>%
  filter(variable == "EOS.pct" & FDR < bal.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_eos_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BE_",.)

genes_eos_be_dn0.3 <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff2 & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif FeNO genes - UP ####
genes_feno_bal_up <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < bal.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_up <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

genes_feno_be_up0.3 <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff2 & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Signif FeNO genes - DOWN ####
genes_feno_bal_dn <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < bal.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

genes_feno_be_dn0.3 <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff2 & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

#### Significant intersects ####
genes_all_bal_up <- intersect(intersect(genes_sbp_bal_up,genes_eos_bal_up),
                              genes_feno_bal_up)
genes_all_be_up <- intersect(intersect(genes_sbp_be_up,genes_eos_be_up),
                             genes_feno_be_up)

genes_all_bal_dn <- intersect(intersect(genes_sbp_bal_dn,genes_eos_bal_dn),
                              genes_feno_bal_dn)
genes_all_be_dn <- intersect(intersect(genes_sbp_be_dn,genes_eos_be_dn),
                             genes_feno_be_dn)

genes_all_be_up0.3 <- intersect(intersect(genes_sbp_be_up,genes_eos_be_up0.3),
                             genes_feno_be_up0.3)
genes_all_be_dn0.3 <- intersect(intersect(genes_sbp_be_dn,genes_eos_be_dn0.3),
                             genes_feno_be_dn0.3)

#### Save ####
save(genes_all_bal_up, genes_all_be_up,
     genes_all_bal_dn, genes_all_be_dn,
     genes_all_be_up0.3, genes_all_be_dn0.3,
     file="publication/signif_genes.RData")
