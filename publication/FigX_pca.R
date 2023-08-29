library(tidyverse)
# library(BIGpicture)
library(patchwork)

#### Data ####
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

dat.BAL.abund.norm.voom$targets <- dat.BAL.abund.norm.voom$targets %>% 
  mutate(visit = fct_recode(visit, "Pre SBP-Ag"="V4", "Post SBP-Ag"="V5"))
dat.BE.abund.norm.voom$targets <- dat.BE.abund.norm.voom$targets %>% 
  mutate(visit = fct_recode(visit, "Pre SBP-Ag"="V4", "Post SBP-Ag"="V5"))

dat.combined <- list()

dat.combined$targets <- bind_rows(dat.BAL.abund.norm.voom$targets,
                                  dat.BE.abund.norm.voom$targets)
dat.combined$E <- inner_join(
  rownames_to_column(as.data.frame(dat.BAL.abund.norm.voom$E)),
  rownames_to_column(as.data.frame(dat.BE.abund.norm.voom$E))) %>% 
  column_to_rownames()

class(dat.combined) <- class(dat.BAL.abund.norm.voom)

#### PCA ####
PCA <- prcomp(t(dat.combined$E), scale. = TRUE, center = TRUE)

PC1.label <- paste("PC1", " (", 
                   summary(PCA)$importance[2, 1] * 100, "%)", sep = "")
PC2.label <- paste("PC2", " (", 
                   summary(PCA)$importance[2, 2] * 100, "%)", sep = "")

pca.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>% 
  left_join(dat.combined$targets)

p1 <- ggplot(pca.dat, aes(PC1,PC2, color = visit, shape = group)) + 
  geom_point(size = 3) + 
  theme_classic() + 
  labs(x = PC1.label, y = PC2.label, shape="", color="") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#74add1", "#FF8679")) +
  theme(legend.position = "bottom", legend.direction = "vertical")
# p1

#### Separate samples ####
##BAL
dat.BAL <- dat.combined
dat.BAL$targets <- dat.BAL$targets %>% filter(group == "BAL")
dat.BAL$E <- dat.BAL$E %>% select(all_of(dat.BAL$targets$libID))

PCA2 <- prcomp(t(dat.BAL$E), scale. = TRUE, center = TRUE)

PC1.label2 <- paste("PC1", " (", 
                   summary(PCA2)$importance[2, 1] * 100, "%)", sep = "")
PC2.label2 <- paste("PC2", " (", 
                   summary(PCA2)$importance[2, 2] * 100, "%)", sep = "")

pca.dat2 <- as.data.frame(PCA2$x) %>% 
  rownames_to_column("libID") %>% 
  left_join(dat.combined$targets)

p2 <- ggplot(pca.dat2, aes(PC1,PC2, color = visit)) + 
  geom_point(size = 3) + 
  theme_classic() + 
  labs(x = PC1.label2, y = PC2.label2, color="", title="BAL") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#74add1", "#FF8679")) +
  theme(legend.position = "none")
# p2

##BE
dat.BE <- dat.combined
dat.BE$targets <- dat.BE$targets %>% filter(group == "BE")
dat.BE$E <- dat.BE$E %>% select(all_of(dat.BE$targets$libID))

PCA3 <- prcomp(t(dat.BE$E), scale. = TRUE, center = TRUE)

PC1.label3 <- paste("PC1", " (", 
                    summary(PCA3)$importance[2, 1] * 100, "%)", sep = "")
PC2.label3 <- paste("PC2", " (", 
                    summary(PCA3)$importance[2, 2] * 100, "%)", sep = "")

pca.dat3 <- as.data.frame(PCA3$x) %>% 
  rownames_to_column("libID") %>% 
  left_join(dat.combined$targets)

p3 <- ggplot(pca.dat3, aes(PC1,PC2, color = visit)) + 
  geom_point(size = 3, shape="triangle") + 
  theme_classic() + 
  labs(x = PC1.label3, y = PC2.label3, color="", title="BE") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#74add1", "#FF8679")) +
  theme(legend.position = "none")
# p3

#### Save ####
plot_all <- p1+p2+p3+ plot_annotation(tag_levels = "A")
plot_all
ggsave(filename = "publication/FigX_pca.png", plot_all,
       width = 8, height = 4)
