library(tidyverse)
library(patchwork)
library(ggrepel)

CDH26_lab <- "CDH26 expression in BE\nlog2 CPM"
cyto_lab <- "IL9 abundance in BAL\nlog10 pg/ml"
eos_lab <- "Eosinophils in BAL (%)"

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
cyto_OI <- c("IL5","IL13","TARC","CCL26","GMCSF","IL9")
multip.V5.be.select <- multip.V5.be[cyto_OI,] %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("sampID")

# Combine
dat <- full_join(multip.V5.be.select, be.cdh26.v5) %>% 
  inner_join(dat.BE$targets %>% 
               filter(visit=="Post") %>% 
               select(donorID, EOS.pct))

#### ACME ####
propmed <- med_result %>% 
  filter(name!="Total" & outcome == 'EOS.pct' &
           name =="Prop_Med") %>% 
  select(name, estimate, cyto) %>% 
  mutate(estimate = case_when(name=="Prop_Med"~estimate*100,
                              TRUE~estimate)) %>% 
  pivot_wider(values_from = estimate)

p1 <- med_result %>% 
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
       color="Percent\nmediated") +
  theme(legend.position = "bottom")
# p1

#### EOS ~ CDH26 ####
pval1 <- med_lm_eos %>%  
  filter(model=="main") %>% 
  distinct(p.value) %>% 
  pull(p.value) %>% unique() %>% 
  signif(., digits=2)

p2 <- ggplot(dat, aes(x=BE_CDH26, y=EOS.pct, color=BE_CDH26)) +
  geom_point() +
  geom_smooth(method="lm", color='black', se=FALSE) +
  theme_classic() +
  labs(x=CDH26_lab,
       y=eos_lab, title="Direct effect") +
  facet_wrap(~paste0("P = ", pval1)) +
  viridis::scale_color_viridis(option = "magma") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
# p2

#### IL9 ~ CDH26 ####
pval2 <- med_lm_eos %>%
  filter(model=="med_out") %>%
  select(cyto, p.value) %>%
  mutate(p.value = signif(p.value, digits=2))

p3 <- dat %>% 
  pivot_longer(IL5:IL9, names_to = "cyto") %>% 
  filter(cyto=="IL9") %>% 
  left_join(pval2) %>%
  mutate(facet.lab = paste0("P = ",p.value)) %>%
  
  ggplot(aes(x=BE_CDH26, y=value, color=BE_CDH26)) +
  geom_point() +
  geom_smooth(method="lm", color='black', se=FALSE) +
  theme_classic() +
  labs(x=CDH26_lab,
       y=cyto_lab,
       title="Mediation effect") +
  facet_wrap(~facet.lab)+
  viridis::scale_color_viridis(option = "magma") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
# p3

#### mediation ####
pval3 <- med_lm_eos %>%
  filter(model=="med" & term=="cyto_value") %>%
  select(cyto, p.value) %>%
  mutate(p.value = signif(p.value, digits=2))

p4 <- dat %>% 
  pivot_longer(IL5:IL9, names_to = "cyto") %>% 
  filter(cyto=="IL9") %>% 
  left_join(pval3) %>%
  mutate(facet.lab = paste0("P = ",p.value)) %>%
  
  ggplot(aes(x=value, y=EOS.pct)) +
  geom_point(aes(color=BE_CDH26)) +
  geom_smooth(method="lm", color='black', se=FALSE) +
  theme_classic() +
  labs(x=cyto_lab,
       y=eos_lab,
       color=CDH26_lab,
       title="Mediator effect") +
  facet_wrap(~facet.lab, nrow=4, scales="free_x") +
  # ylim(c(-30,90)) +
  viridis::scale_color_viridis(option = "magma") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(hjust = 0.5), 
        legend.title.position = "bottom")
# p4

#### Change in beta ####
# #Summarise change in estimate
# 
# dat_beta <- med_lm_eos %>%  
#   # filter(cyto %in% cyto_OI) %>% 
#   filter(model=="med" & term == "BE_CDH26") %>% 
#   #add non-mediator row
#   bind_rows(med_lm_eos %>%  
#               # filter(cyto %in% cyto_OI) %>% 
#               filter(model=="main" & term == "BE_CDH26") %>% 
#               distinct(outcome,model, term, estimate, p.value) %>% 
#               mutate(cyto="None")) %>% 
#   left_join(propmed) %>% 
#   mutate(cyto = factor(cyto),
#          cyto = fct_relevel(cyto, "None", after=0L))
# 
# p5 <- dat_beta %>% 
#   ggplot(aes(x=reorder(cyto,-estimate), 
#              y=estimate, fill=Prop_Med)) +
#   geom_col() +
#   theme_bw() +
#   viridis::scale_fill_viridis() +
#   geom_hline(yintercept = max(dat_beta$estimate),
#              lty="dashed") +
#   labs(x="Mediating cytokine",
#        y="CDH26 expression effects\non eosinophil % (B)",
#        fill="Percent\nmediated")  +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1))

#### cytokine correlation ####
cor.dat <- multip.V5.be.select %>% 
  column_to_rownames("sampID") %>% 
  as.matrix() 
cor.result <- Hmisc::rcorr(cor.dat, type = "pearson")

#Arrange as in barplot
cyto_ord <-rev(c("IL9","IL13","IL5","TARC","CCL26","GMCSF"))

cor.arrange <- as.data.frame(cor.result$r) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "r")

cor.arrange2 <- as.data.frame(cor.result$P) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "p") %>% 
  full_join(cor.arrange) %>% 
  #order Y
  mutate(name = factor(name, levels=cyto_ord),
         X_variable = factor(X_variable, 
                             levels=cyto_ord)) %>% 
  arrange(name, X_variable) %>% 
  rowwise() %>%
  mutate(pair = sort(c(X_variable, name)) %>% 
           paste(collapse = ",")) %>%
  #Remove lower triangle
  group_by(pair) %>%
  distinct(pair, .keep_all = T) %>%
  filter(X_variable != name) %>%
  ungroup() %>%
  #signif label
  mutate(label = ifelse(p < 0.01, "*",""),
         FDR = p.adjust(p),
         label2 = ifelse(FDR < 0.01, "*",""))

p6 <- as.data.frame(cor.arrange2) %>% 
  ggplot(aes(x=X_variable, y=name, fill=r)) +
  geom_tile() +
  # geom_text(aes(label = label)) +
  labs(x="", fill="Pearson R") +
  coord_fixed()+
  # coord_flip() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_gradient(low="white",high="darkred",
                       limits=c(0.5,1), breaks=c(0.5,0.75,1.0))
# p6

#### Save ####
# lo <-"
# IABD
# EFCH
# "
# p_all <- p1+p2+p3+p4+p5+p6 + plot_spacer()+plot_spacer()+
#   plot_layout(design = lo, widths = 1, heights = 1)+
#   plot_annotation(tag_levels = list(c("A","E","F","G","C","D","","B")))
lo <-"
ACE
BDF
"
p_all <- p1+p6+p2+p3+p4 +plot_spacer()+
  plot_layout(design = lo, widths = 1, heights = 1)+
  plot_annotation(tag_levels = list(c("A","B","C"," "," "," ","")))

# p_all

ggsave("publication/Fig4_CDH26_c.pdf", p_all,
       width=9.1, height=7)


