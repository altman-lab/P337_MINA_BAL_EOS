library(tidyverse)
library(Hmisc)
library(corrplot)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

samp.overlap <- intersect(unique(dat.BAL.abund.norm.voom$targets$donorID),
                          unique(dat.BE.abund.norm.voom$targets$donorID))

#format
targets.bal <- dat.BAL.abund.norm.voom$targets %>% 
  filter(donorID %in% samp.overlap) %>% 
  rownames_to_column() %>% 
  distinct(donorID, visit, EOS.pct, NEUT.pct, LYM.pct, MONO.pct,Epi.pct) %>% 
  rename(Eosinophil=EOS.pct, Neutrophil=NEUT.pct, Monocyte=MONO.pct,
         Epithelial=Epi.pct, Lymphocyte=LYM.pct)

targets.bal.v4 <- targets.bal %>% 
  filter(visit == "V4") %>% 
  select(-visit) %>% 
  column_to_rownames("donorID") %>% 
  as.matrix()
targets.bal.v5 <- targets.bal %>% 
  filter(visit == "V5") %>% 
  select(-visit) %>% 
  column_to_rownames("donorID") %>% 
  as.matrix()

#### Corr ####
corr.v4 <- rcorr(targets.bal.v4, type = "pearson")
corr.v5 <- rcorr(targets.bal.v5, type = "pearson")

#### Heatmap ####
png("publication/FigX_cell.corr.png", width=8, height=5, units = "in", res = 300)
par(mfrow=c(1,2))
corrplot(corr.v4$r, p.mat = corr.v4$P, method = 'color', 
         diag = FALSE, type = 'lower', col = rev(COL2("RdBu")),
         sig.level = c(0.01, 0.1), pch.cex = 1.5,
         insig = 'label_sig',
         tl.col = "black", title = "Pre SBP-Ag", mar=c(0,0,1,0))
corrplot(corr.v5$r, p.mat = corr.v5$P, method = 'color', 
         diag = FALSE, type = 'lower', col = rev(COL2("RdBu")),
         sig.level = c(0.01, 0.1), pch.cex = 1.5,
         insig = 'label_sig',
         tl.col = "black", title = "Post SBP-Ag", mar=c(0,0,1,0))
dev.off()
