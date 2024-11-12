library(tidyverse)
library(patchwork)
library(ggrepel)

CDH26_lab <- "CDH26 expression in BE (log2 CPM)"
cyto_lab <- "Cytokine abundance in BAL (log10 pg/ml)"
eos_lab <- "Eosinophil % in BAL"

#### Data ####
#Gene expression
attach("data_clean/P337_BE_data.RData")
dat.BE$targets <- dat.BE$targets %>% 
  mutate(visit = case_when(visit == "V4" ~ "Pre",
                           visit == "V5" ~ "Post"),
         visit = factor(visit, levels = c("Pre","Post")))

# Cytokine protein abundance
multip <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  #Correct error
  rename(FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into = c("name","visit"), remove = FALSE) %>% 
  mutate(visit = case_when(visit == "V4" ~ "Pre",
                           visit == "V5" ~ "Post")) %>% 
  mutate(sampID = paste(ptID,visit,sep="_")) %>% 
  select(-ptID,-visit) %>% 
  mutate(value=log10(value)) %>% 
  pivot_wider(names_from = name) %>% 
  arrange(sampID) %>% 
  filter(grepl("_Post", sampID)) %>% 
  column_to_rownames("sampID") %>% 
  t() %>% as.data.frame()

#Mediation
attach("results/mediation.RData")
med_lm_eos <- med_lm_result %>% 
  filter(term!="(Intercept)") %>% 
  select(outcome, cyto, model, term, estimate, p.value) %>%
  filter(outcome=="EOS.pct") 

#### CDH26 ####
cdh26ens <- dat.BE$genes %>% 
  filter(hgnc_symbol=="CDH26") %>% 
  pull(geneName)

be.cdh26 <- as.data.frame(dat.BE$E) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% cdh26ens) %>% 
  pivot_longer(-rowname, names_to = "libID", values_to="BE_CDH26") %>% 
  left_join(dat.BE$targets) %>% 
  mutate(sampID = paste(donorID,visit,sep="_")) %>% 
  select(sampID, donorID, visit, BE_CDH26) %>% 
  filter(visit=="Post") %>% 
  arrange(sampID)

#Overlapping samples
be_overlap <- intersect(be.cdh26$sampID, colnames(multip))
multip.V5.be <- multip %>% 
  select(all_of(be_overlap))
be.cdh26.v5 <- be.cdh26 %>% 
  filter(sampID %in% be_overlap)

# Select cytokines of interest
cyto_OI <- c("IL5","IL13","TARC","IL9")
multip.V5.be.select <- multip.V5.be[cyto_OI,] %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("sampID")

# Combine
dat <- full_join(multip.V5.be.select, be.cdh26.v5) %>% 
  inner_join(dat.BE$targets %>% 
              filter(visit=="Post") %>% 
              select(donorID, EOS.pct))

#### EOS ~ CDH26 ####
pval1 <- med_lm_eos %>%  
  filter(model=="main") %>% 
  distinct(p.value) %>% 
  pull(p.value) %>% unique() %>% 
  signif(., digits=2)

# est1 <- med_lm_eos %>%  
#   filter(model=="main") %>% 
#   distinct(estimate) %>% 
#   pull(estimate) %>% unique() %>% 
#   signif(., digits=3)

p1 <- ggplot(dat, aes(x=BE_CDH26, y=EOS.pct)) +
  geom_point(color="grey30") +
  geom_smooth(method="lm", color='black') +
  theme_classic() +
  labs(x=gsub("in","\nin",CDH26_lab),
       y=eos_lab) +
  facet_wrap(~paste0("CDH26\nP = ", pval1)) +
  ylim(c(-30,90))
# p1

#### Cyto ~ CDH26 ####
pval2 <- med_lm_eos %>%
  filter(model=="med_out") %>%
  select(cyto, p.value) %>%
  mutate(p.value = signif(p.value, digits=2))

# est2 <- med_lm_eos %>%
#   filter(model=="med_out") %>%
#   select(cyto, estimate) %>%
#   mutate(estimate = signif(estimate, digits=2))

p2 <- dat %>% 
  pivot_longer(IL5:IL9, names_to = "cyto") %>% 
  left_join(pval2) %>%
  # left_join(est2) %>%
  mutate(facet.lab = paste0(cyto,#"\nestimate = ",estimate,
                            "\nP = ",p.value)) %>%
  
  ggplot(aes(x=BE_CDH26, y=value)) +
  geom_point(color="grey30") +
  geom_smooth(method="lm", color='black') +
  theme_classic() +
  labs(x=CDH26_lab,
       y=cyto_lab) +
  facet_wrap(~facet.lab, nrow=2)
# p2

#### mediation ####
pval3 <- med_lm_eos %>%
  filter(model=="med" & term=="cyto_value") %>%
  select(cyto, p.value) %>%
  mutate(p.value = signif(p.value, digits=2))

# est3 <- med_lm_eos %>%
#   filter(model=="med" & term=="cyto_value") %>%
#   select(cyto, estimate) %>%
#   mutate(estimate = signif(estimate, digits=2))

p3 <- dat %>% 
  pivot_longer(IL5:IL9, names_to = "cyto") %>% 
  left_join(pval3) %>%
  # left_join(est3) %>%
  mutate(facet.lab = paste0(cyto,#"\nestimate = ",estimate,
                            "\nP = ",p.value)) %>%
  
  ggplot(aes(x=value, y=EOS.pct)) +
  geom_point(aes(color=BE_CDH26)) +
  geom_smooth(method="lm", color='black') +
  theme_classic() +
  labs(x=gsub("in ","\nin ",cyto_lab),
       y=eos_lab,
       color=gsub("in","\nin",CDH26_lab)) +
  facet_wrap(~facet.lab, nrow=4, scales="free_x") +
  ylim(c(-30,90)) +
  viridis::scale_color_viridis() +
  theme(legend.position = "bottom")
  
# p3

#### combo effect: EOS ~ CDH26 + cyto ####
#Summarise change in estimate and ACME

med_lm_eos %>%  
  filter(cyto %in% cyto_OI) %>% 
  filter(model=="med" & term == "BE_CDH26") %>% 
  arrange(cyto)


#### ACME ####
propmed <- med_result %>% 
  filter(name!="Total" & outcome == 'EOS.pct' &
           name =="Prop_Med") %>% 
  select(name, estimate, cyto) %>% 
  mutate(estimate = case_when(name=="Prop_Med"~estimate*100,
                              TRUE~estimate)) %>% 
  pivot_wider(values_from = estimate)

p4 <- med_result %>% 
  filter(name!="Total" & outcome == 'EOS.pct' &
           name =="ACME") %>% 
  select(cyto, name, estimate, p) %>% 
  full_join(propmed) %>% 
  
  # mutate(p=ifelse(p==0, 2E-16, p)) %>% 
  ggplot(aes(x=estimate, y=-log10( p ), color=Prop_Med)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), lty="dashed") +
  geom_text_repel(aes(label = cyto), 
                  show.legend = FALSE, max.overlaps = 100,
                  color="black")+
  theme_bw() +
  viridis::scale_color_viridis() +
  labs(x = "Average causal mediation effect\n(ACME)",
       y="-log10( p )",
       color="Percent mediated") 
# p4


#### Save ####
lo <-"
AA
BC
DD
DD
EF
EF
EF
EF
"
p_all <- p4+p1+plot_spacer()+p2+p3+plot_spacer()+
  plot_layout(design = lo)+
  plot_annotation(tag_levels = "A")
# p_all

ggsave("publication/Fig4_CDH26.pdf", p_all,
       width=5, height=12)


