library(tidyverse)
library(BIGpicture)
select <- dplyr::select

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("results/model_fitting.RData")
bal.cutoff <- 0.1
be.cutoff <- 0.3

## Signif EOS genes - UP
genes_eos_bal_up <- lm_eos_BAL$lm %>%
  filter(variable == "EOS.pct" & FDR < bal.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_eos_be_up <- lm_eos_BE$lm %>%
  filter(variable == "EOS.pct" & FDR < be.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

## Signif EOS genes - DOWN
genes_eos_bal_dn <- lm_eos_BAL$lm %>%
  filter(variable == "EOS.pct" & FDR < bal.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BAL_",.)

genes_eos_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate < 0) %>%
  pull(gene) %>% unique() %>% paste0("BE_",.)

## Signif FeNO genes - UP
genes_feno_bal_up <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < bal.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_up <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate > 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

## Signif FeNO genes - DOWN
genes_feno_bal_dn <- lm_feno_BAL$lm %>%
  filter(variable == "FeNO" & FDR < bal.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BAL_",.) 

genes_feno_be_dn <- lm_feno_BE$lm %>%
  filter(variable == "FeNO" & FDR < be.cutoff & estimate < 0) %>% 
  pull(gene) %>% unique() %>% paste0("BE_",.)

## Significant intersects
genes_both_bal_up <- intersect(genes_eos_bal_up, 
                               c(genes_feno_bal_up,genes_feno_bal_dn))
genes_both_be_up <- intersect(genes_eos_be_up, 
                              c(genes_feno_be_up,genes_feno_be_dn))

genes_both_bal_dn <- intersect(genes_eos_bal_dn, 
                               c(genes_feno_bal_up,genes_feno_bal_dn))
genes_both_be_dn <- intersect(genes_eos_be_dn, 
                              c(genes_feno_be_up,genes_feno_be_dn))

#### Enrichment for coloring ####
all_genes <- c(genes_both_bal_up,genes_both_be_up,
               genes_both_bal_dn,genes_both_be_dn)

enrich <- data.frame(gene=all_genes) %>% 
  mutate(pathway = c(rep("BAL increase eosinophils and FeNO",
                         length(genes_both_bal_up)),
                     rep("BE increase eosinophils and FeNO",
                         length(genes_both_be_up)),
                     rep("BAL decrease eosinophils and FeNO",
                         length(genes_both_bal_dn)),
                     rep("BE decrease eosinophils and FeNO",
                         length(genes_both_be_dn)))) %>% 
  separate(gene, c("group","geneName"), sep = "_") %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  group_by(pathway) %>% 
  summarise(genes=list(hgnc_symbol)) %>% 
  mutate(group="BAL", FDR=0, `k/K`=1, group_in_pathway=10)
  
#### string ####
all_genes_hgnc <- enrich %>% 
  unnest(genes) %>% 
  pull(genes) %>% unique()

map <- map_string(genes = unique(all_genes_hgnc),
                  score_threshold = 400)
  
p1 <- plot_string(map, enrichment = enrich, layout = "fr",
            colors = c("#0571b0","#ca0020","#92c5de","#f4a582"))

ggsave(p1, filename="publication/FigX_string2.png",
       width=15, height=10)

#### enrichment ####
library(SEARchways)

enrich.dat <- enrich %>% 
  unnest(genes) %>% 
  select(pathway, genes) %>% 
  dplyr::rename(group=pathway, gene=genes)
  
H <- BIGprofiler(gene_df = enrich.dat, category = "H", ID = "SYMBOL")

p4 <- plot_enrich(H, fdr_cutoff = 0.1)
p4

ggsave(p4, filename="publication/FigX_enrich.png", 
       width=14, height=3)

p2 <- plot_string(map, enrichment = H, layout = "fr")
# p2


####
GO <- BIGprofiler(gene_df = enrich.dat, category = "C5",
                  subcategory = "GO:BP",
                  ID = "SYMBOL")

# p3 <- plot_string(map, enrichment = GO, layout = "fr",
#                   fdr_cutoff = 0.1)
# p3

# plot_enrich(GO, fdr_cutoff = 0.1)
library(rrvgo)
library(GO.db)

# add GO ID
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`) %>% 
  mutate(GO_path_format = gsub("-|,|\\/"," ",tolower(GO_path)),
         GO_path_format = gsub("  "," ", GO_path_format))

enrich.anno <- GO %>% 
  filter(FDR <= 0.2 & group_in_pathway >1) %>%
  mutate(GO_path_format = gsub("GOBP_","",pathway),
         GO_path_format = gsub("_"," ",tolower(GO_path_format))) %>%
  left_join(goterms)

rrvgo.ls <- list()
for(g in unique(enrich.anno$group)){
  temp <- enrich.anno %>% filter(group==g)
  
  if(nrow(temp)>1){
    simMatrix <- calculateSimMatrix(temp$GOID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")
    
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    
    rrvgo.ls[[g]] <- reducedTerms
  }
}

rrvgo.ls7 <- rrvgo.ls
treemapPlot(rrvgo.ls7[[1]], title = names(rrvgo.ls7)[1])
treemapPlot(rrvgo.ls7[[2]], title = names(rrvgo.ls7)[2])
treemapPlot(rrvgo.ls7[[3]], title = names(rrvgo.ls7)[3])
# treemapPlot(rrvgo.ls7[[4]], title = names(rrvgo.ls7)[4])
