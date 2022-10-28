
########################################
##                                    ##
##          Install Packages          ##
##                                    ##
########################################

install.packages('Seurat')
install.packages("spatstat")

library(Seurat)

# Needs python
library(reticulate)
py_config()
# create a new environment 
conda_create("r-reticulate")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

install.packages("UpSetR")
install.packages('VennDiagram')

### You will need this to update R ###
#install.packages("installr")
#install.packages("stringr")
#library(installr)
#install.packages("dplyr")
#install.packages("cowplot")
#updateR()

########################################
##                                    ##
##           Load Packages            ##
##                                    ##
########################################

library(Seurat)
library(cowplot)
library(dplyr)
library(sctransform)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(scales)
library(viridis)
library(viridisLite)
library(UpSetR)
library(VennDiagram)

########################################
##                                    ##
##      Code for Making Figures       ##
##                                    ##
########################################
load("XXX/PBMC_Immune.combined_scT_every_filt_20191017.Robj")
setwd("XXX/RFigs")

# Example
load("C:/Users/adkim/Documents/NKCellTests/PBMC_Immune.combined_scT_every_filt_20191017.Robj")
setwd("C:/Users/adkim/Documents/NKCellTests")


# This code here just reorganizes the cell types for aesthetic reasons. 
# I keep it all here so you can organize it easily however you want, but you'll need to read the labels

immune.combined$type <- factor(immune.combined$type, levels = c("Healthy","Alcohol"))

# preserved is just the label for the original Seurat clusters that I then renamed somewhat generically eg Cytotoxic_T-cell1
Idents(immune.combined) <- "preserved"
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("CD4_T-cell1","CD4_T-cell2","CD4_T-cell3",
                                                                      "Cytotoxic_T-cell1", "Cytotoxic_T-cell2","Cytotoxic_T-cell3","Cytotoxic_T-cell4",
                                                                      "CD14_Monocyte1","CD14_Monocyte2","CD14_Monocyte3","CD16_Monocyte",
                                                                      "B-cell","NK-cell",
                                                                      "pDC","Megakaryocyte",
                                                                      "Unassigned"))
immune.combined$celltype <- factor(immune.combined$celltype, levels = c("CD4_T-cell1_Healthy","CD4_T-cell2_Healthy","CD4_T-cell3_Healthy",
                                                                        "Cytotoxic_T-cell1_Healthy", "Cytotoxic_T-cell2_Healthy",
                                                                        "Cytotoxic_T-cell3_Healthy","Cytotoxic_T-cell4_Healthy",
                                                                        "CD14_Monocyte1_Healthy","CD14_Monocyte2_Healthy",
                                                                        "CD14_Monocyte3_Healthy","CD16_Monocyte_Healthy",
                                                                        "B-cell_Healthy","NK-cell_Healthy",
                                                                        "pDC_Healthy","Megakaryocyte_Healthy",
                                                                        "Unassigned_Healthy",
                                                                        "CD4_T-cell1_Alcohol","CD4_T-cell2_Alcohol","CD4_T-cell3_Alcohol",
                                                                        "Cytotoxic_T-cell1_Alcohol", "Cytotoxic_T-cell2_Alcohol",
                                                                        "Cytotoxic_T-cell3_Alcohol","Cytotoxic_T-cell4_Alcohol",
                                                                        "CD14_Monocyte1_Alcohol","CD14_Monocyte2_Alcohol",
                                                                        "CD14_Monocyte3_Alcohol","CD16_Monocyte_Alcohol",
                                                                        "B-cell_Alcohol","NK-cell_Alcohol",
                                                                        "pDC_Alcohol","Megakaryocyte_Alcohol",
                                                                        "Unassigned_Alcohol",
                                                                        "CD4_T-cell1_NA","CD4_T-cell2_NA","CD4_T-cell3_NA",
                                                                        "Cytotoxic_T-cell1_NA", "Cytotoxic_T-cell2_NA",
                                                                        "Cytotoxic_T-cell3_NA","Cytotoxic_T-cell4_NA",
                                                                        "CD14_Monocyte1_NA","CD14_Monocyte2_NA",
                                                                        "CD14_Monocyte3_NA","CD16_Monocyte_NA",
                                                                        "B-cell_NA","NK-cell_NA",
                                                                        "pDC_NA","Megakaryocyte_NA",
                                                                        "Unassigned_NA"))
immune.combined$celltype.stim <- factor(immune.combined$celltype.stim, levels = c("CD4_T-cell1_Healthy_Basal","CD4_T-cell2_Healthy_Basal","CD4_T-cell3_Healthy_Basal",
                                                                        "Cytotoxic_T-cell1_Healthy_Basal", "Cytotoxic_T-cell2_Healthy_Basal",
                                                                        "Cytotoxic_T-cell3_Healthy_Basal","Cytotoxic_T-cell4_Healthy_Basal",
                                                                        "CD14_Monocyte1_Healthy_Basal","CD14_Monocyte2_Healthy_Basal",
                                                                        "CD14_Monocyte3_Healthy_Basal","CD16_Monocyte_Healthy_Basal",
                                                                        "B-cell_Healthy_Basal","NK-cell_Healthy_Basal",
                                                                        "pDC_Healthy_Basal","Megakaryocyte_Healthy_Basal",
                                                                        "Unassigned_Healthy_Basal",
                                                                        "CD4_T-cell1_Alcohol_Basal","CD4_T-cell2_Alcohol_Basal","CD4_T-cell3_Alcohol_Basal",
                                                                        "Cytotoxic_T-cell1_Alcohol_Basal", "Cytotoxic_T-cell2_Alcohol_Basal",
                                                                        "Cytotoxic_T-cell3_Alcohol_Basal","Cytotoxic_T-cell4_Alcohol_Basal",
                                                                        "CD14_Monocyte1_Alcohol_Basal","CD14_Monocyte2_Alcohol_Basal",
                                                                        "CD14_Monocyte3_Alcohol_Basal","CD16_Monocyte_Alcohol_Basal",
                                                                        "B-cell_Alcohol_Basal","NK-cell_Alcohol_Basal",
                                                                        "pDC_Alcohol_Basal","Megakaryocyte_Alcohol_Basal",
                                                                        "Unassigned_Alcohol_Basal"))


#Only for DimPLot
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("CD14_Monocyte1","CD14_Monocyte2","CD14_Monocyte3","CD16_Monocyte",
                                                                      "CD4_T-cell1","CD4_T-cell2","CD4_T-cell3",
                                                                      "Cytotoxic_T-cell1", "Cytotoxic_T-cell2","Cytotoxic_T-cell3","Cytotoxic_T-cell4",
                                                                      "B-cell","NK-cell",
                                                                      "pDC","Megakaryocyte",
                                                                      "Unassigned"))
DimPlot(immune.combined, label = TRUE, label.size = 5)

# This was the original UMAP plot for all the data that was published in Hep Communications: 
# DOI: 10.1002/hep4.1563
pdf("Clustering_OldFig.pdf", height = 5, width = 7, useDingbats=FALSE)
p <- DimPlot(immune.combined, label = TRUE, label.size = 3)
p
dev.off()


########################################
##                                    ##
##       Baseline Only Figures        ##
##           Figure 1D                ##
########################################
Idents(immune.combined) <- "celltype.stim"

DefaultAssay(immune.combined) <- "SCT"

pdf("Violin_CD8sub_KLRB1_KLRG1_noLPS.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("KLRB1", "KLRG1"), split.by = "type", 
                 group.by = "preserved", split.plot = TRUE,pt.size = 0, combine = FALSE, 
                 idents = c('NK-cell_Healthy_Basal','NK-cell_Alcohol_Basal',
                            'Cytotoxic_T-cell1_Healthy_Basal','Cytotoxic_T-cell1_Alcohol_Basal',
                            'Cytotoxic_T-cell2_Healthy_Basal','Cytotoxic_T-cell2_Alcohol_Basal',
                            'Cytotoxic_T-cell3_Healthy_Basal','Cytotoxic_T-cell3_Alcohol_Basal',
                            'Cytotoxic_T-cell4_Healthy_Basal','Cytotoxic_T-cell4_Alcohol_Basal'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'grey'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("Violin_CD8sub_KLRC1_KLRC2_noLPS.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("KLRC1", "KLRC2"), split.by = "type", 
                 group.by = "preserved", split.plot = TRUE,pt.size = 0, combine = FALSE, 
                 idents = c('NK-cell_Healthy_Basal','NK-cell_Alcohol_Basal',
                            'Cytotoxic_T-cell1_Healthy_Basal','Cytotoxic_T-cell1_Alcohol_Basal',
                            'Cytotoxic_T-cell2_Healthy_Basal','Cytotoxic_T-cell2_Alcohol_Basal',
                            'Cytotoxic_T-cell3_Healthy_Basal','Cytotoxic_T-cell3_Alcohol_Basal',
                            'Cytotoxic_T-cell4_Healthy_Basal','Cytotoxic_T-cell4_Alcohol_Basal'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'grey'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("Violin_CD8sub_GZMA_GZMB_noLPS.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("GZMA", "GZMB"), split.by = "type", 
                 group.by = "preserved", split.plot = TRUE,pt.size = 0, combine = FALSE, 
                 idents = c('NK-cell_Healthy_Basal','NK-cell_Alcohol_Basal',
                            'Cytotoxic_T-cell1_Healthy_Basal','Cytotoxic_T-cell1_Alcohol_Basal',
                            'Cytotoxic_T-cell2_Healthy_Basal','Cytotoxic_T-cell2_Alcohol_Basal',
                            'Cytotoxic_T-cell3_Healthy_Basal','Cytotoxic_T-cell3_Alcohol_Basal',
                            'Cytotoxic_T-cell4_Healthy_Basal','Cytotoxic_T-cell4_Alcohol_Basal'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'grey'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("Violin_CD8sub_PRF1_GNLY_noLPS.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(immune.combined, features = c("PRF1", "GNLY"), split.by = "type", 
                 group.by = "preserved", split.plot = TRUE,pt.size = 0, combine = FALSE, 
                 idents = c('NK-cell_Healthy_Basal','NK-cell_Alcohol_Basal',
                            'Cytotoxic_T-cell1_Healthy_Basal','Cytotoxic_T-cell1_Alcohol_Basal',
                            'Cytotoxic_T-cell2_Healthy_Basal','Cytotoxic_T-cell2_Alcohol_Basal',
                            'Cytotoxic_T-cell3_Healthy_Basal','Cytotoxic_T-cell3_Alcohol_Basal',
                            'Cytotoxic_T-cell4_Healthy_Basal','Cytotoxic_T-cell4_Alcohol_Basal'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('black', 'grey'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()



########################################
##                                    ##
##         Gene Correlation           ##
##                                    ##
########################################

#BIGscale

#Install and open packages
install.packages('devtools')
library(devtools)
install.packages('ellipsis')

# Make sure these get downloaded - sometimes it fails and you must try again
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("BioQC")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# This is a fussy download, try a few times, and don't update things it asks for, usually unnecessary
devtools::install_github("iaconogi/bigSCale2")

library("bigSCale", lib.loc="C:/Users/adkim/AppData/Local/R/win-library/4.2")

# Convert Seurat data to SingeCellExperiment Type
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html

install.packages("BiocManager")
BiocManager::install("scater")
install.packages("textshape")

# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')

library("bigSCale", lib.loc="C:/Users/adkim/AppData/Local/R/win-library/4.2")
library(devtools)
library(scater)
library(bigSCale)
library(loomR)
library(Seurat)
library(textshape)
library(ggplot2)
library(reshape2)
library(igraph)
library(ggplot2)
library(reshape2)
library(cowplot)

load("XXX/PBMC_Immune.combined_scT_every_filt_20191017.Robj")
setwd("XXX/RFigs")

# Example
load("C:/Users/adkim/Documents/NKCellTests/PBMC_Immune.combined_scT_every_filt_20191017.Robj")
setwd("C:/Users/adkim/Documents/NKCellTests")

# List of Compiled DE Genes

targets <- read.csv("XXX/All_genes_v2.csv", header = FALSE)

targets <- read.csv("C:/Users/adkim/Documents/NKCellTests/All_genes_v2.csv", header = FALSE)
NKCtargets <-targets$V1 
rm(targets)

###################
# HC Basal Figure #
###################
HC_B_NK_c <- subset(immune.combined, idents = c("NK-cell_Healthy_Basal"))
HC_B_NK_c <- FindVariableFeatures(HC_B_NK_c, selection.method = "vst", nfeatures = 2000)

DefaultAssay(HC_B_NK_c) <- "SCT"
HC_B.data <- as.data.frame(as.matrix(GetAssayData(HC_B_NK_c, slot = "counts")))
gene.names <- rownames(HC_B_NK_c)
HC_B.results=compute.network(expr.data = HC_B.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
HC_B.results$graph
DT::datatable(HC_B.results$centrality)

HCBmatrix<-as.data.frame(as.numeric(HC_B.results$correlations))

# You are my density
allHCB<-melt(HCBmatrix)
#HCB <- plot(density(allHCB$value))

matrixtest <- HCBmatrix
matrixall_rank_all <- allHCB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed4 <-cluster_matrix(matrix_subbed3, dim='both', method='complete')
matrix_subbed4$gene <- rownames(matrix_subbed4)

#d <- dist(matrix_subbed3, method = "euclidean")
#hc1 <- hclust(d, method = "complete" )
#plot(hc1, cex = 0.6, hang = -1)
#write.table(hc1$order, file="Healthy_Basal_HierClust.txt")

NKC_all <- melt(matrix_subbed4)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed4$gene))

NKC_rank_all <- NKC_all

d <- density(NKC_rank_all$value)
plot(d)

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

write.table(NKC_rank_all, file="NKCorr_HC_Basal.txt")

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

# This creates Figure 1E
pdf("HCB_AllDE.pdf", height = 14, width = 14, useDingbats=FALSE)
p
dev.off()

#################
# HC LPS Figure #
#################
HC_L_NK_c <- subset(immune.combined, idents = c("NK-cell_Healthy_LPS"))
HC_L_NK_c <- FindVariableFeatures(HC_L_NK_c, selection.method = "vst", nfeatures = 2000)

DefaultAssay(HC_L_NK_c) <- "SCT"
HC_L.data <- as.data.frame(as.matrix(GetAssayData(HC_L_NK_c, slot = "counts")))
gene.names <- rownames(HC_L_NK_c)
HC_L.results=compute.network(expr.data = HC_L.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
HC_L.results$graph
DT::datatable(HC_L.results$centrality)

HCLmatrix<-as.data.frame(as.numeric(HC_L.results$correlations))

# You are my density
allHCL<-melt(HCLmatrix)
#HCL <- plot(density(allHCL$value))

matrixtest <- HCLmatrix
matrixall_rank_all <- allHCL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed4 <-cluster_matrix(matrix_subbed3, dim='both', method='complete')
matrix_subbed4$gene <- rownames(matrix_subbed4)

#d <- dist(matrix_subbed3, method = "euclidean")
#hc1 <- hclust(d, method = "complete" )
#plot(hc1, cex = 0.6, hang = -1)
#write.table(hc1$order, file="Healthy_Basal_HierClust.txt")

NKC_all <- melt(matrix_subbed4)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed4$gene))

NKC_rank_all <- NKC_all

d <- density(NKC_rank_all$value)
plot(d)

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

write.table(NKC_rank_all, file="NKCorr_HC_LPS.txt")

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

# This creates Figure 2B
pdf("HCL_AllDE.pdf", height = 14, width = 14, useDingbats=FALSE)
p
dev.off()

###################
# AH Basal Figure #
###################
AH_B_NK_c <- subset(immune.combined, idents = c("NK-cell_Alcohol_Basal"))
AH_B_NK_c <- FindVariableFeatures(AH_B_NK_c, selection.method = "vst", nfeatures = 2000)

DefaultAssay(AH_B_NK_c) <- "SCT"
AH_B.data <- as.data.frame(as.matrix(GetAssayData(AH_B_NK_c, slot = "counts")))
gene.names <- rownames(AH_B_NK_c)
AH_B.results=compute.network(expr.data = AH_B.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
AH_B.results$graph
DT::datatable(AH_B.results$centrality)

AHBmatrix<-as.data.frame(as.numeric(AH_B.results$correlations))

# You are my density
allAHB<-melt(AHBmatrix)
#AHB <- plot(density(allAHB$value))

matrixtest <- AHBmatrix
matrixall_rank_all <- allAHB

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed4 <-cluster_matrix(matrix_subbed3, dim='both', method='complete')
matrix_subbed4$gene <- rownames(matrix_subbed4)

#d <- dist(matrix_subbed3, method = "euclidean")
#hc1 <- hclust(d, method = "complete" )
#plot(hc1, cex = 0.6, hang = -1)
#write.table(hc1$order, file="Healthy_Basal_HierClust.txt")

NKC_all <- melt(matrix_subbed4)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed4$gene))

NKC_rank_all <- NKC_all

d <- density(NKC_rank_all$value)
plot(d)

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

write.table(NKC_rank_all, file="NKCorr_AH_Basal.txt")

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

# This creates Figure 1F
pdf("AHB_AllDE.pdf", height = 14, width = 14, useDingbats=FALSE)
p
dev.off()

#################
# AH LPS Figure #
#################
AH_L_NK_c <- subset(immune.combined, idents = c("NK-cell_Alcohol_LPS"))
AH_L_NK_c <- FindVariableFeatures(AH_L_NK_c, selection.method = "vst", nfeatures = 2000)

DefaultAssay(AH_L_NK_c) <- "SCT"
AH_L.data <- as.data.frame(as.matrix(GetAssayData(AH_L_NK_c, slot = "counts")))
gene.names <- rownames(AH_L_NK_c)
AH_L.results=compute.network(expr.data = AH_L.data,
                             gene.names = gene.names
                             #clustering = 'direct'
)
AH_L.results$graph
DT::datatable(AH_L.results$centrality)

AHLmatrix<-as.data.frame(as.numeric(AH_L.results$correlations))

# You are my density
allAHL<-melt(AHLmatrix)
#AHL <- plot(density(allAHL$value))

matrixtest <- AHLmatrix
matrixall_rank_all <- allAHL

keep <- as.character(NKCtargets)
matrix_subbed <- matrixtest[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))
matrix_subbed3 <- matrix_subbed2[keep,]
matrix_subbed3 <- na.omit(matrix_subbed3)
matrix_subbed4 <-cluster_matrix(matrix_subbed3, dim='both', method='complete')
matrix_subbed4$gene <- rownames(matrix_subbed4)

#d <- dist(matrix_subbed3, method = "euclidean")
#hc1 <- hclust(d, method = "complete" )
#plot(hc1, cex = 0.6, hang = -1)
#write.table(hc1$order, file="Healthy_Basal_HierClust.txt")

NKC_all <- melt(matrix_subbed4)
NKC_all$gene <- factor(NKC_all$gene, levels = c(matrix_subbed4$gene))

NKC_rank_all <- NKC_all

d <- density(NKC_rank_all$value)
plot(d)

NKC_rank_all$rank <-0
NKC_rank_all$rank[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]<- NKC_rank_all$value[which(NKC_rank_all$value>quantile(matrixall_rank_all$value, 0.975))]
NKC_rank_all$rank[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]<- NKC_rank_all$value[which(NKC_rank_all$value<quantile(matrixall_rank_all$value, 0.025))]

write.table(NKC_rank_all, file="NKCorr_AH_LPS.txt")

p <- ggplot(NKC_rank_all, aes(gene, variable)) + 
  geom_tile(aes(fill = rank), colour = "white") + 
  theme(axis.text.x = element_text(size = 10 ,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  theme(legend.position = "none")

# This creates Figure 2C
pdf("AHL_AllDE.pdf", height = 14, width = 14, useDingbats=FALSE)
p
dev.off()




########################################
##                                    ##
##            Upset Test              ##
##                                    ##
########################################

setwd("XXX/RFigs")
NK_CD8_LPS_UP_Reorder <- read.csv("###/NK_CD8_LPS_UP_Reorder.csv")

# Example
setwd("C:/Users/adkim/Documents/NKCellTests")
NK_CD8_LPS_UP_Reorder <- read.csv("C:/Users/adkim/Documents/NKCellTests/NK_CD8_LPS_UP_Reorder.csv")

# This will create Figure 2A
list <- colnames(NK_CD8_LPS_UP_Reorder)
upset(NK_CD8_LPS_UP_Reorder, 
      sets = list[2:7], 
      #group.by = "sets", # This will reorder the bars, turn off for published version
      cutoff = 5,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Genes Per Cell Type", 
      #empty.intersections = "on",
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 1.5), # First number is y axis label, Last number is number above bar
      order.by = "freq")


#############################################
###                                       ###
###               Venn Diagram            ###
###                                       ###
#############################################

# FYI - This script always results in +1 added to the middle (the total intersection)
# I have never solved this problem out of laziness, 
# but it's probably due to an extra line at the end of each column, 
# giving an intersection of null for all columns

##############################
# This will create Figure 1B #
##############################
NK_CD8_Basal_UP <- read.csv("XXX/NK_CD8_Basal_UP.csv")
setwd("XXX/RFigs")

# Example
NK_CD8_Basal_UP <- read.csv("C:/Users/adkim/Documents/NKCellTests/NK_CD8_Basal_UP.csv")
setwd("C:/Users/adkim/Documents/NKCellTests")

Z<-calculate.overlap(list(NK=NK_CD8_Basal_UP[,1], CD8_1=NK_CD8_Basal_UP[,2], 
                          CD8_2=NK_CD8_Basal_UP[,3], CD8_3=NK_CD8_Basal_UP[,4], CD8_4=NK_CD8_Basal_UP[,5]))
a<-venn.diagram(x=list(NK=NK_CD8_Basal_UP[,1], CD8_1=NK_CD8_Basal_UP[,2], 
                       CD8_2=NK_CD8_Basal_UP[,3], CD8_2=NK_CD8_Basal_UP[,4], CD8_2=NK_CD8_Basal_UP[,5]), 
                filename = "NK_Basal_Up.tiff")


##############################
# This will create Figure 1C #
##############################
NK_CD8_Basal_UP <- read.csv("XXX/NK_CD8_Basal_DOWN.csv")
setwd("XXX/RFigs")

# Example
NK_CD8_Basal_DOWN <- read.csv("C:/Users/adkim/Documents/NKCellTests/NK_CD8_Basal_DOWN.csv")
setwd("C:/Users/adkim/Documents/NKCellTests")

Z<-calculate.overlap(list(NK=NK_CD8_Basal_UP[,1], CD8_1=NK_CD8_Basal_UP[,2], 
                          CD8_2=NK_CD8_Basal_UP[,3] ))
a<-venn.diagram(x=list(NK=NK_CD8_Basal_DOWN[,1], CD8_1=NK_CD8_Basal_DOWN[,2], 
                       CD8_2=NK_CD8_Basal_DOWN[,3]), 
                filename = "NK_Basal_Down.tiff")


