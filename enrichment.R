#### Enrichment ####
library(tidyverse)
library(ComplexHeatmap)
library(SEARchways)
source("publication/get_hm_cluster.R")

load("publication/heatmap_dat.RData")

#### UP ####
#Cell type
c7_up <- BIGprofiler(gene_df = hm_clust_up, ID = "SYMBOL",
                  category = "C7", subcategory = "IMMUNESIGDB")

c7_up_signif <- c7_up %>% 
  filter(grepl("GSE3982_", pathway)) %>% 
  mutate(FDR = p.adjust(pval, method="BH")) %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>% 
  select(group, pathway, group_in_pathway, `k/K`, pval, FDR) %>% 
  separate(pathway, into=c("study","rest"), sep="_", extra="merge", 
           remove = FALSE) %>% 
  separate(rest, into=c("cell1","cell2"), sep="_VS_", extra="merge") %>% 
  separate(cell2, into=c("cell2","direction"), sep="_(?!.*_)", extra="merge") %>% 
  mutate(cell_assoc = case_when(direction == "UP"~ cell1,
                                direction == "DN"~ cell2)) %>% 
  filter(cell_assoc != "CTRL")

c7_up_signif %>% count(group) 

c7_up_signif %>% 
  count(group, cell_assoc) %>% 
  arrange(group)%>% 
  rename(`GSE3982 terms`=n)

#GO
c5_up <- BIGprofiler(gene_df = hm_clust_up, ID = "SYMBOL",
                  category = "C5", subcategory = "GO:BP")

c5_up_signif <- c5_up %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>% 
  select(group, pathway, group_in_pathway, `k/K`, pval, FDR) 

c5_up_signif %>% count(group) 

# c5_up_signif2 <- c5_up %>% 
#   filter(FDR < 0.2 & group_in_pathway > 1) %>% 
#   select(group, pathway, pathway_GOID, group_in_pathway, FDR, genes) 
# 
# c5_up_signif2 %>% count(group)

#### DOWN ####
#Cell type
c7_dn <- BIGprofiler(gene_df = hm_clust_dn, ID = "SYMBOL",
                     category = "C7", subcategory = "IMMUNESIGDB")

c7_dn_signif <- c7_dn %>% 
  filter(grepl("GSE3982_", pathway)) %>% 
  mutate(FDR = p.adjust(pval, method="BH")) %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>% 
  select(group, pathway, group_in_pathway, `k/K`, pval, FDR) %>% 
  separate(pathway, into=c("study","rest"), sep="_", extra="merge", 
           remove = FALSE) %>% 
  separate(rest, into=c("cell1","cell2"), sep="_VS_", extra="merge") %>% 
  separate(cell2, into=c("cell2","direction"), sep="_(?!.*_)", extra="merge") %>% 
  mutate(cell_assoc = case_when(direction == "UP"~ cell1,
                                direction == "DN"~ cell2)) %>% 
  filter(cell_assoc != "CTRL")

c7_dn_signif %>% count(group) 

c7_dn_signif %>% 
  count(group, cell_assoc) %>% 
  arrange(group)%>% 
  rename(`GSE3982 terms`=n)

#GO
c5_dn <- BIGprofiler(gene_df = hm_clust_dn, ID = "SYMBOL",
                     category = "C5", subcategory = "GO:BP")

c5_dn_signif <- c5_dn %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>% 
  select(group, pathway, group_in_pathway, `k/K`, pval, FDR) 

c5_dn_signif %>% count(group) 

# c5_dn_signif2 <- c5_dn %>% 
#   filter(FDR < 0.2 & group_in_pathway > 1) %>% 
#   select(group, pathway, pathway_GOID, group_in_pathway, FDR, genes) 
# 
# c5_dn_signif2 %>% count(group)

#### Save ####
save(c7_up_signif, c5_up_signif, c7_dn_signif, c5_dn_signif,
     file="results/DEG_enrichment.RData")

#### RRVGO ####
library(rrvgo)
plot.ls <- list()

for(clust in sort(unique(c5.signif2$group))){
  print(clust)
  if(clust == "cluster4"){
    #### FDR < 0.2 ####
    temp2 <- c5.signif2 %>% filter(group==clust)
    
    simMatrix2 <- calculateSimMatrix(temp2$pathway_GOID,
                                     orgdb="org.Hs.eg.db",
                                     ont="BP",
                                     method="Rel")
    
    scores2 <- setNames(-log10(temp2$FDR), temp2$pathway_GOID)
    reducedTerms2 <- reduceSimMatrix(simMatrix2,
                                     scores2,
                                     threshold=0.8,
                                     orgdb="org.Hs.eg.db")
    plot.ls[[paste(clust,"FDR<0.2")]] <- treemapPlot(reducedTerms2, 
                                                     title=paste(clust,"FDR<0.2"))
  } else {
    #### FDR < 0.1 ####
    temp1 <- c5.signif1 %>% filter(group==clust)
    
    simMatrix1 <- calculateSimMatrix(temp1$pathway_GOID,
                                     orgdb="org.Hs.eg.db",
                                     ont="BP",
                                     method="Rel")
    
    scores1 <- setNames(-log10(temp1$FDR), temp1$pathway_GOID)
    reducedTerms1 <- reduceSimMatrix(simMatrix1,
                                     scores1,
                                     threshold=0.8,
                                     orgdb="org.Hs.eg.db")
    plot.ls[[paste(clust,"FDR<0.1")]] <- treemapPlot(reducedTerms1, 
                                                     title=paste(clust,"FDR<0.1"))
  }
}
