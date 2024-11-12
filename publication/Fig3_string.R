library(tidyverse)
library(BIGpicture)
library(patchwork)
select <- dplyr::select

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("results/gene_model_fitting.RData")
load("publication/signif_genes.RData")

#### Coloring ####
all_genes <- c(genes_all_bal_up,genes_all_be_up,
               genes_all_bal_dn,genes_all_be_dn)

enrich <- data.frame(gene=all_genes) %>% 
  mutate(pathway = c(rep("BAL positive",
                         length(genes_all_bal_up)),
                     rep("BE positive",
                         length(genes_all_be_up)),
                     rep("BAL negative",
                         length(genes_all_bal_dn)),
                     rep("BE negative",
                         length(genes_all_be_dn)))) %>% 
  separate(gene, c("group","geneName"), sep = "_") %>% 
  left_join(dat.BAL$genes) %>% 
  group_by(pathway) %>% 
  summarise(genes=list(hgnc_symbol)) %>% 
  mutate(group="BAL", FDR=0, `k/K`=1, group_in_pathway=10)
  
#### Map to string ####
all_genes_hgnc <- enrich %>% 
  unnest(genes) %>% 
  pull(genes) %>% unique()
all_genes_unique <- sort(unique(all_genes_hgnc))

map <- map_string(genes = all_genes_unique,
                  score_threshold = 150)

#### string : everything ####
plot_string(map, enrichment = enrich, layout = "fr",
                  #edge_min = 1,
            colors = c("#92c5de","#f4a582","#0571b0","#ca0020"))

#### Separate cluster from other ####
genes_remove <- c("CEP68","DCAF4","HSD17B8","TMEM74B","TMEM150B","RUSC1")

map2 <- map_string(genes = all_genes_unique[!all_genes_unique %in% genes_remove],
                   score_threshold = 150)
map3 <- map_string(genes = genes_remove,
                   score_threshold = 150)

p1 <- plot_string(map2, enrichment = enrich, layout = "fr",
                  colors = c("#92c5de","#f4a582","#0571b0","#ca0020"))
p1

#### small cluster and indivs ####
p2 <- plot_string(map3, enrichment = enrich, layout = "grid",
                  colors = c("#92c5de","#0571b0","#ca0020"))
p2

#### Save ####
plot_all <- p1/p2
plot_all

ggsave(plot_all, filename="publication/Fig3_string150.pdf",
       width=8, height=14)
