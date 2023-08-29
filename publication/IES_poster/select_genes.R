library(tidyverse)
library(readxl)
library(patchwork)
# library(ggpubr)
library(Hmisc)

#### Targeted genegenes ####
genes.OI <- c("MUC5AC", "MUC5B", "FOXA3","GALNT6", "TFF1", "TFF2", "TFF3","LYZ")

#### Data ####
# BAL data
load("data_clean/P337_BAL_data.RData")
# BE data
load("data_clean/P337_BE_data.RData")

#### Combine ####
temp <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  #gene expression
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  mutate(group = "BAL") %>% 
  #sample metadata
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              select(libID,donorID,visit)) %>% 
  #gene names
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  filter(hgnc_symbol %in% genes.OI)

dat <- as.data.frame(dat.BE.abund.norm.voom$E) %>% 
  #gene expression
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  mutate(group = "BE") %>% 
  #sample metadata
  left_join(dat.BE.abund.norm.voom$targets %>% 
              select(libID,donorID,visit)) %>% 
  #gene names
  left_join(dat.BE.abund.norm.voom$genes) %>% 
  filter(hgnc_symbol %in% genes.OI) %>% 
  bind_rows(temp) 

dat_format <- dat %>% 
  select(-libID) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(diff = V5-V4,
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) %>% 
  pivot_longer(V4:V5, names_to = "visit") %>% 
  left_join(dat.BE.abund.norm.voom$targets %>% 
              select(libID,donorID,visit)) %>% 
  #fix labels
  mutate(visit=recode_factor(visit, "V4"="Pre", "V5"="Post"))%>% 
  drop_na(diff)

#### Pre vs post plots ####
p1 <- dat_format %>% 
  filter(group=="BAL") %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial provocation with allergen (SBP-Ag)",
       y="Gene log2 CPM", color="Post - Pre\nchange") +
  facet_wrap(~hgnc_symbol, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679"))
# p1

p2 <- dat_format %>% 
  filter(group=="BE") %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial provocation with allergen (SBP-Ag)",
       y="Gene log2 CPM", color="Post - Pre\nchange") +
  facet_wrap(~hgnc_symbol, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679"))
# p2

#### delta corr to L vA insula ####
#fMRI data
neuro <- read_excel(sheet="T1",
                    "../P337_MINA_BAL_public/data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add visit variable
  mutate(visit="V4")

#Time point 2 = visit 5
neuro <- read_excel(sheet="T2",
                    "../P337_MINA_BAL_public/data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add visit variable
  mutate(visit="V5") %>% 
  #Combine with other visit
  full_join(neuro) %>% 
  #Format idnum to match RNAseq data
  mutate(idnum = paste("MA",idnum, sep="")) %>% 
  rename(donorID=idnum)

neuro.delta <- neuro %>% 
  filter(donorID != "MA1012" & neuro != "LSI" & neuro != "BDI") %>% 
  #Calculate delta
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  #remove fxnal metrics
  filter(!grepl("3way",neuro) & !grepl("LPR", neuro)) %>% 
  #fix name
  mutate(neuro = gsub("amgy","amyg",neuro)) %>% 
  #wide format
  dplyr::select(-V4,-V5) %>% 
  pivot_wider(names_from = neuro, values_from = delta) 

#genes
plot_ls <- list()
plot_ls2 <- list()

for (g in sort(unique(dat_format$hgnc_symbol))){
  print(g)
  bal_temp <- dat_format %>%  filter(group=="BAL")
  be_temp <- dat_format %>%  filter(group=="BE")
  
  if(g %in% unique(bal_temp$hgnc_symbol)){
    gene_temp <- bal_temp %>% 
      filter(group=="BAL") %>% 
      filter(hgnc_symbol==g) %>% 
      select(-libID) %>% 
      pivot_wider(names_from = visit) %>% 
      mutate(delta=Post-Pre)  %>% 
      select(-Pre,-Post) %>% 
      left_join(neuro.delta %>% select(donorID, L_vA_insula),
                by="donorID") %>% 
      drop_na(L_vA_insula, delta)
    
    corr0 <- rcorr(x=gene_temp$L_vA_insula,
                   y=gene_temp$delta, type="pearson") 
    
    title0 <- paste(paste("R", signif(corr0$r[1,2], digits=2), sep=" = "),
                    paste("P", signif(corr0$P[1,2], digits=2), sep=" = "), sep=", ")
    
    if(g==sort(unique(dat_format$hgnc_symbol))[1]){
      ylab <-"Post - pre L vA insula"
    } else{
      ylab <- ""
    }
    
    plot_ls[[g]] <- gene_temp %>% 
      ggplot(aes(y=L_vA_insula, x=delta)) +
      geom_point() +
      geom_smooth(method="lm",se=FALSE, color="#FF8679", formula = 'y ~ x') +
      labs(y=ylab,
           x=paste("Post - pre", g, "\nlog2 CPM"),
           title=title0) +
      theme_classic()+
      theme(plot.title  = element_text(size=10))
  }
  
  ## BE plot
  if(g %in% unique(be_temp$hgnc_symbol)){
    gene_temp2 <- be_temp %>% 
      filter(hgnc_symbol==g) %>% 
      select(-libID) %>% 
      pivot_wider(names_from = visit) %>% 
      mutate(delta=Post-Pre)  %>% 
      select(-Pre,-Post) %>% 
      left_join(neuro.delta %>% select(donorID, L_vA_insula),
                by="donorID") %>% 
      drop_na(L_vA_insula, delta)
    
    corr02 <- rcorr(x=gene_temp2$L_vA_insula,
                    y=gene_temp2$delta, type="pearson") 
    
    title02 <- paste(paste("R", signif(corr02$r[1,2], digits=2), sep=" = "),
                     paste("P", signif(corr02$P[1,2], digits=2), sep=" = "), sep=", ")
    
    plot_ls2[[g]] <- gene_temp2 %>% 
      ggplot(aes(y=L_vA_insula, x=delta)) +
      geom_point() +
      geom_smooth(method="lm",se=FALSE, color="#FF8679", formula = 'y ~ x') +
      labs(y=ylab,
           x=paste("Post - pre", g, "\nlog2 CPM"),
           title=title0) +
      theme_classic()+
      theme(plot.title  = element_text(size=10))
  }
}


corr_plots1 <- wrap_plots(plot_ls, nrow = 1, tag_level="new")
corr_plots2 <- wrap_plots(plot_ls2, nrow = 1, tag_level="new")

#### Save ####

plot_BAL <- 
  p1 / corr_plots1 + plot_annotation(title="BAL")
plot_BE <- 
  p2 / corr_plots2 + plot_annotation(title="BE")

ggsave(plot_BAL, filename="publication/FigX.BAL.more.genes.png",
       width=11.5, height=7)
ggsave(plot_BE, filename="publication/FigX.BE.more.genes.png",
       width=14, height=7)
