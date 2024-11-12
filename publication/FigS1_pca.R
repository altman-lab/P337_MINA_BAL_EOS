library(tidyverse)
# library(BIGpicture)
library(patchwork)
library(edgeR)

#### Data ####
#final data
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BE_data.RData")

dat.combined <- list()

dat.combined$targets <- bind_rows(dat.BAL$targets,
                                  dat.BE$targets)
dat.combined$E <- inner_join(
  rownames_to_column(as.data.frame(dat.BAL$E)),
  rownames_to_column(as.data.frame(dat.BE$E))) %>% 
  column_to_rownames()

class(dat.combined) <- class(dat.BAL)

#Pre batch filtering
load("data_raw/counts.pc.filter.RData")

#### Batch ####
#Calculate PCA for all data.
PCA.all <- counts.pc.filter %>% 
  column_to_rownames("geneName") %>% 
  #Convert to log counts per million
  cpm(., log=TRUE) %>% 
  t() %>% 
  #Calc PCA
  prcomp()

PC1.label <- paste("PC1 (", summary(PCA.all)$importance[2,1]*100, "%)", sep="")
PC2.label <-paste("PC2 (", summary(PCA.all)$importance[2,2]*100, "%)", sep="")

# Extract PC values
PCA.all.dat <- as.data.frame(PCA.all$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs for plotting
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  full_join(meta.filter, by="libID") %>% 
  mutate(`SBP-Ag` = recode(visit, "V4"="Pre","V5"="Post"))

PCA.batch <- PCA.all.dat %>% 
  
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color=paste(group, batch), shape=`SBP-Ag`),
             size=3) + 
  #Beautify
  theme_classic() +
  theme(legend.title = element_blank()) +
  labs(x=PC1.label, y=PC2.label, 
       title="Un-normalized log2CPM") +
  coord_fixed(ratio=1) + 
  scale_color_manual(values=c("#882255","#88CCEE","#DDCC77","#44AA99"))

#list duplicate samples remaining in data
dups <- PCA.all.dat %>% 
  mutate(dupID = paste(group, donorID, visit, sep="_")) %>% 
  select(dupID) %>% unlist(use.names = FALSE)
dups <- dups[duplicated(dups)]

PCA.dups <- PCA.all.dat %>% 
  mutate(dupID = paste(group, donorID, visit, sep="_")) %>% 
  mutate(col.group = ifelse(dupID %in% dups, "duplicate", "none")) %>% 
  arrange(desc(col.group)) %>% 
  
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color=col.group, shape=`SBP-Ag`),
             size=3) + 
  geom_path(aes(group=dupID)) +
  #Beautify
  theme_classic() +
  theme(legend.title = element_blank()) +
  labs(x=PC1.label, y=PC2.label, 
       title="Un-normalized log2CPM") +
  coord_fixed(ratio=1) +
  scale_color_manual(values=c("#CC6677","#117733"))
# PCA.dups

PCA.batch/PCA.dups + 
  plot_layout(guides = "collect")

#### PCA ####
PCA <- prcomp(t(dat.combined$E), scale. = TRUE, center = TRUE)

PC1.label1 <- paste("PC1", " (", 
                   summary(PCA)$importance[2, 1] * 100, "%)", sep = "")
PC2.label1 <- paste("PC2", " (", 
                   summary(PCA)$importance[2, 2] * 100, "%)", sep = "")

pca.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>% 
  left_join(dat.combined$targets)

p1 <- ggplot(pca.dat, aes(PC1,PC2, color = visit, shape = group)) + 
  geom_point(size = 3) + 
  theme_classic() + 
  labs(x = PC1.label1, y = PC2.label1, shape="", color="",
       title="Normalized log2CPM") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#74add1", "#FF8679"))

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
  labs(x = PC1.label2, y = PC2.label2, color="", 
       title="BAL normalized log2CPM") + 
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
  labs(x = PC1.label3, y = PC2.label3, color="", 
       title="BE normalized log2CPM") + 
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#74add1", "#FF8679")) +
  theme(legend.position = "none")
# p3

#### Save ####
plot_all <- (PCA.batch+PCA.dups)/(p1+p2+p3)+ 
  plot_annotation(tag_levels = "A")
plot_all

ggsave(filename = "publication/FigS1_pca.png", plot_all,
       width = 12, height = 7)
ggsave(filename = "publication/FigS1_pca.pdf", plot_all,
       width = 12, height = 7)
