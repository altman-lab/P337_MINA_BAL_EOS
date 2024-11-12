library(tidyverse)
library(patchwork)
library(ggpubr)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")
attach("results/gene_model_fitting.RData")

#gene expression
dat.BAL.rename <- as.data.frame(dat.BAL$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BAL$targets %>% 
              select(libID, donorID, visit, EOS.pct, FeNO.PreBro_V5)) %>% 
  left_join(dat.BAL$genes, by=c("gene"="geneName")) %>% 
  select(hgnc_symbol, visit, donorID, value, EOS.pct, FeNO.PreBro_V5) %>% 
  mutate(group="BAL")

dat.BE.rename <- as.data.frame(dat.BE$E) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.BE$targets %>% 
              select(libID, donorID, visit, EOS.pct, FeNO.PreBro_V5)) %>% 
  left_join(dat.BE$genes, by=c("gene"="geneName")) %>% 
  select(hgnc_symbol, visit, donorID, value, EOS.pct, FeNO.PreBro_V5) %>% 
  mutate(group="BE")

#Combine bal and be expression
dat.all.rename <- dat.BE.rename %>% 
  bind_rows(dat.BAL.rename) %>% 
  arrange(donorID) %>% 
  mutate(visit = fct_recode(visit, "Pre"="V4", "Post"="V5"))

#LM results
lm.bal <- bind_rows(lm_visit_BAL$lm, lm_eos_BAL$lm, lm_feno_BAL$lm) %>% 
  mutate(group = "BAL") %>% 
  left_join(dat.BAL$genes, by=c("gene"="geneName")) %>% 
  select(group, hgnc_symbol, variable, FDR)

lm.be <- bind_rows(lm_visit_BE$lm, lm_eos_BE$lm, lm_feno_BE$lm) %>% 
  mutate(group = "BE") %>% 
  left_join(dat.BE$genes, by=c("gene"="geneName")) %>% 
  select(group, hgnc_symbol, variable, FDR)

lm.all <- bind_rows(lm.bal, lm.be) %>% 
  mutate(FDR = paste0(group, "\n","FDR = ", signif(FDR, digits=2))) %>% 
  filter(variable != "(Intercept)") %>% 
  mutate(variable = paste0(variable,"_FDR")) 

#### Plots ####
geneOI <- c("SERPINB2")
plot.ls <- list()

for(g in geneOI){
  dat.temp <- dat.all.rename %>% 
    filter(hgnc_symbol==g)
  
  lm.temp <- lm.all %>% 
    filter(hgnc_symbol==g) %>% 
    pivot_wider(names_from = variable, values_from = FDR) %>% 
    full_join(dat.temp, by=c("group","hgnc_symbol"))
  
  #### SBP plot ####
  p1 <- ggplot(lm.temp, aes(x=visit, y=value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0) +
    theme_classic() +
    facet_wrap(~visit_FDR) +
    labs(x="Segmented bronchial\nprovaction with allergen\n(SBP-Ag)",
         y="Normalized log2 expression")
  
  #### EOS plot ####
  p2 <- lm.temp %>% 
    filter(visit=="Post") %>% 
    ggplot(aes(x=EOS.pct, y=value)) +
    geom_point() +
    geom_smooth(formula="y~x", method = 'lm', color="grey", se=FALSE) +
    theme_classic() +
    facet_wrap(~EOS.pct_FDR, scales="free", ncol=1) +
    labs(x="eosinophils (%)",
         y="Normalized log2 expression",
         title="Post SBP-Ag")
  
  #### FEV1 plot ####
  p3 <- lm.temp %>% 
    filter(visit=="Post") %>% 
    ggplot(aes(x=FeNO.PreBro_V5, y=value)) +
    geom_point() +
    geom_smooth(formula="y~x", method = 'lm', color="grey", se=FALSE) +
    theme_classic() +
    facet_wrap(~FeNO_FDR, scales = "free", ncol=1) +
    labs(x="FeNO",
         y="Normalized log2 expression",
         title = "Post SBP-Ag")
  
  #### BAL BE corr ####
  p4 <- lm.temp %>% 
    filter(visit=="Post") %>% 
    select(group,value, donorID) %>% 
    pivot_wider(names_from = group) %>% 
    ggplot(aes(x=BAL, y=BE)) +
    geom_point() +
    geom_smooth(formula="y~x", method = 'lm', color="grey", se=FALSE) +
    stat_cor(method = "pearson") +
    theme_classic() +
    labs(x="BAL Post SBP-Ag\nNormalized log2 expression",
         y="BE Post SBP-Ag\nNormalized log2 expression") +
    coord_fixed()
  
  lo <- "
  ABC
  ABC
  ABC
  ADD
  ADD
  "
  p_all <- p1 + p2 + p3 + p4 + plot_layout(design = lo)+
    plot_annotation(tag_levels = "A",
                    title=g)
  
  ggsave(p_all,
         filename = paste0("publication/FigX_",g,".png"), 
         width=9, height=7)
  plot.ls[[g]] <- p_all
}


#### Save ####

