library(tidyverse)
library(limma)
library(ggpubr)
library(patchwork)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

samp.overlap <- intersect(unique(dat.BAL.abund.norm.voom$targets$donorID),
                          unique(dat.BE.abund.norm.voom$targets$donorID))

#format
targets.all <- dat.BAL.abund.norm.voom$targets %>% 
  filter(donorID %in% samp.overlap) %>% 
  distinct(donorID, visit, EOS.pct, NEUT.pct, 
           FEV1.pctPP.preAlbuterol_V4, FEV1.pctPP.preAlbuterol_V5, 
           FeNO.PreBro_V4, FeNO.PreBro_V5) %>% 
  pivot_longer(-c(donorID, visit, EOS.pct, NEUT.pct), values_to = "value") %>% 
  separate(name, into=c("name","visit2"), sep="_") %>% 
  filter(visit==visit2) %>% 
  pivot_wider() %>% 
  rename(FEV1=FEV1.pctPP.preAlbuterol, FeNO=FeNO.PreBro)

#### Cells ####
targets.all.format <- targets.all %>% 
  distinct(donorID, visit, EOS.pct, NEUT.pct) %>% 
  pivot_longer(-c(donorID, visit), names_to = "variable") %>% 
  mutate(fdr = recode(variable, 
                      "EOS.pct"="eosinophils\nFDR = 9.1E-5",
                      "NEUT.pct"="neutrophils\nFDR = 2.0E-2"),
         visit = fct_recode(visit, "Pre"="V4", "Post"="V5")) %>% 
  group_by(donorID, variable) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) 


p1 <- targets.all.format %>% 
  filter(variable %in% c("EOS.pct","NEUT.pct")) %>% 
  ggplot() +
  aes(x=visit, y= value) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  #Beautify
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(linewidth=0.5),
        legend.position = "bottom",
        legend.direction = "vertical") +
  labs(y="Cells (%)", x="Segmented bronchial\nprovocation with allergen\n(SBP-Ag)",
       color="Post - Pre change") +
  facet_wrap(~fdr, scales="free_x", ncol=1) +
  scale_color_manual(values = c("#D55E00","#0072B2"))

# p1

#### Corr to FeNO ####

dat.eos <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, EOS.pct) %>% 
  mutate(name="EOS.pct") %>% 
  rename(value=EOS.pct)

dat.feno <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, FeNO.PreBro_V4:FeNO.PreBro_V5) %>% 
  pivot_longer(-donorID) %>%
  separate(name, into = c("name","visit"), sep="_") %>% 
  distinct()

p2 <- full_join(dat.eos,dat.feno) %>% 
  pivot_wider() %>%
  mutate(visit = fct_recode(visit,"Pre SPB-Ag"="V4","Post SBP-Ag"="V5")) %>% 

  ggplot(aes(x=EOS.pct,y=FeNO.PreBro)) +
  geom_point(color="grey") +
  theme_classic() +
  stat_cor() +
  geom_smooth(method="lm", se=FALSE, formula = "y~x", color="black") +
  labs(x = "eosinophils (%)", y = "FeNO") +
  facet_wrap(~visit, scales="free", nrow=1) +
  theme(legend.position = "none")
# p2

dat.neut <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, NEUT.pct) %>% 
  mutate(name="NEUT.pct") %>% 
  rename(value=NEUT.pct)

p3 <- full_join(dat.neut,dat.feno) %>% 
  pivot_wider() %>%
  mutate(visit = fct_recode(visit,"Pre SPB-Ag"="V4","Post SBP-Ag"="V5")) %>% 
  
  ggplot(aes(x=NEUT.pct,y=FeNO.PreBro)) +
  geom_point(color="grey") +
  theme_classic() +
  stat_cor() +
  geom_smooth(method="lm", se=FALSE, formula = "y~x", color="black") +
  labs(x = "neutrophils (%)", y = "FeNO") +
  facet_wrap(~visit, scales="free", nrow=1) +
  theme(legend.position = "none")
# p3


#### Corr to FEV ####

dat.eos <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, EOS.pct) %>% 
  mutate(name="EOS.pct") %>% 
  rename(value=EOS.pct)

dat.fev <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, FEV1.pctPP.preAlbuterol_V4:FEV1.pctPP.preAlbuterol_V5) %>% 
  pivot_longer(-donorID) %>%
  separate(name, into = c("name","visit"), sep="_") %>% 
  distinct()

p4 <- full_join(dat.eos,dat.fev) %>% 
  pivot_wider() %>%
  mutate(visit = fct_recode(visit,"Pre SPB-Ag"="V4","Post SBP-Ag"="V5")) %>% 
  
  ggplot(aes(x=EOS.pct,y=FEV1.pctPP.preAlbuterol)) +
  geom_point(color="grey") +
  theme_classic() +
  stat_cor() +
  geom_smooth(method="lm", se=FALSE, formula = "y~x", color="black") +
  labs(x = "eosinophils (%)", y = "FEV1 (%)") +
  facet_wrap(~visit, scales="free", nrow=1) +
  theme(legend.position = "none")
# p4

dat.neut <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, NEUT.pct) %>% 
  mutate(name="NEUT.pct") %>% 
  rename(value=NEUT.pct)

p5 <- full_join(dat.neut,dat.fev) %>% 
  pivot_wider() %>%
  mutate(visit = fct_recode(visit,"Pre SPB-Ag"="V4","Post SBP-Ag"="V5")) %>% 
  
  ggplot(aes(x=NEUT.pct,y=FEV1.pctPP.preAlbuterol)) +
  geom_point(color="grey") +
  theme_classic() +
  stat_cor() +
  geom_smooth(method="lm", se=FALSE, formula = "y~x", color="black") +
  labs(x = "neutrophils (%)", y = "FEV1 (%)") +
  facet_wrap(~visit, scales="free", nrow=1) +
  theme(legend.position = "none")
# p5


#### Save ####
design <- "
ABBB
ACCC
"

plot_all <- p1+p2+p3 + plot_layout(design = design) + 
  plot_annotation(tag_levels = "A")
# plot_all
ggsave(filename = "publication/FigX_cells.feno.png", plot_all,
       width=8.5, height=7)

plot_all2 <- p4/p5 + 
  plot_annotation(tag_levels = "A")
# plot_all2
ggsave(filename = "publication/FigSX_fev1.png", plot_all2,
       width=6, height=6)
