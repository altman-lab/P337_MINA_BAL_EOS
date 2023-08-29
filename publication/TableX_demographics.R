library(tidyverse)
library(limma)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

samp.overlap <- intersect(unique(dat.BAL.abund.norm.voom$targets$donorID),
                          unique(dat.BE.abund.norm.voom$targets$donorID))

meta.raw <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  filter(donorID %in% samp.overlap) 

##Numeric
meta <- meta.raw  %>% 
  mutate(age_yr = age_mo/12) %>% 
  distinct(donorID, age_yr, `FEV1 screening WLAC`, `FEV1 PP`) %>% 
  pivot_longer(-donorID) %>% 
  mutate(visit="none")

meta2 <- meta.raw %>% 
  select(donorID, FEV1.pctPP.preAlbuterol_V4:FeNO.PreBro_V5) %>% 
  distinct() %>% 
  pivot_longer(-donorID) %>% 
  separate(name, into=c("name","visit"), sep="_")

meta3 <- meta.raw %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  distinct() %>% 
  pivot_longer(-c(donorID, visit))

meta.num <- bind_rows(meta, meta2, meta3) %>% 
  arrange(donorID, name)

##Categorical
meta.cat <- meta.raw %>% 
  distinct(donorID, sex, race) %>% 
  pivot_longer(-donorID)

#### Summary stats ###

meta.raw %>% 
  count(group, visit)

meta.num %>% 
  group_by(name, visit) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm=TRUE)) %>% 
  arrange(visit) %>% View()

meta.num %>%
  pivot_wider(names_from = visit) %>% 
  select(-none) %>% 
  drop_na() %>% 
  group_by(name) %>% 
  summarise(ttest = t.test(V4, V5, paired = TRUE)["p.value"]$p.value) %>% 
  ungroup() %>% 
  mutate(group = ifelse(grepl(".pct$",name),"cell","other")) %>% 
  group_by(group) %>% 
  mutate(FDR = p.adjust(ttest, method="fdr"))

meta.cat %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(total = sum(n),
         perc = n/total*100)

