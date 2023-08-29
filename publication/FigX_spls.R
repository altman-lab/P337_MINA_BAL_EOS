library(tidyverse)
library(tidytext)
library(patchwork)
library(ggpubr)

set.seed(3698)
`%notin%` <- Negate(`%in%`)
select <- dplyr::select

#### PLS Data ####
load("results/PLS/PLS_data.RData")
load("results/PLS/SPLS.RData")
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

samp.overlap <- intersect(unique(dat.BAL.abund.norm.voom$targets$donorID),
                          unique(dat.BE.abund.norm.voom$targets$donorID))

#### Loadings ####
loadX <- as.data.frame(spls.tuned$loadings$X) %>% 
  rownames_to_column() %>% 
  mutate(space="X") %>% 
  pivot_longer(comp1, names_to = "comp") %>% 
  mutate(name = case_when(comp=="comp1"~
                            paste("Component 1 = ",
                                  round(spls.tuned$prop_expl_var$X[1]*100,digits=1), "%"))) %>% 
  separate(rowname, into=c("group","geneName"), sep="_", remove = FALSE) %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  mutate(rowname2 = paste(group, hgnc_symbol)) %>% 
  filter(abs(value)>0)

p1 <- loadX %>% 
  ggplot(aes(x=reorder(rowname2, abs(value)), 
             y=value)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  facet_wrap(~name, scales="free", ncol=2) +
  labs(x="Selected gene", y="Loading", fill="") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0) +
  lims(y=c(-0.75,0.75))
# p1

#### Linear regression ####
p2 <- as.data.frame(X) %>% 
  rownames_to_column("donorID") %>% 
  pivot_longer(-donorID) %>% 
  filter(name %in% loadX$rowname) %>% 
  separate(name, into=c("group","geneName")) %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              filter(visit=="V5") %>%
              select(donorID, EOS.pct)) %>% 

  ggplot(aes(x=EOS.pct,y=value)) +
  geom_point(color="grey") +
  theme_classic() +
  stat_cor() +
  geom_smooth(method="lm", se=FALSE, formula = "y~x", color="black") +
  labs(x = "Eosinophils post SBP-Ag (%)", y = "Log2 normalized expression pre SBP-Ag") +
  facet_wrap(~hgnc_symbol, scales="free", ncol=2) +
  theme(legend.position = "none")
p2

#### Save ####
plot_all <- p1/p2 + plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1,3))
# plot_all

ggsave("publication/FigX_spls.png", plot_all,
       width=5, height=6)
