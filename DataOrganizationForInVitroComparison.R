library(Seurat)
library(tidyverse)

#### HUMAN HYPO - ALL CELLS #### FOR PLOT 1
KaZhouAll <- readRDS("/Users/hannahglover/Dropbox/Columbia/Brian Hypothalamus/ScienceAdvances_2023_Revision2/CLEAN_SEURAT/KaZhouAll.rds")
HumanHypoAllCells <- as.data.frame(rowMeans(KaZhouAll@assays$RNA@counts))
colnames(HumanHypoAllCells) <- "HumanHypo_AllCells_AllTimepoints"
row.names(HumanHypoAllCells) <- row.names(KaZhouAll@assays$RNA) #36116 genes

for(x in unique(KaZhouAll@meta.data$Timepoint)){
  Idents(KaZhouAll) <- "Timepoint"
  KaZhouAll_ByTimepoint <- subset(KaZhouAll, idents = x)
  
  HumanHypoTimepoint <- as.data.frame(rowMeans(KaZhouAll_ByTimepoint@assays$RNA@counts))
  colnames(HumanHypoTimepoint) <- paste0("HYPO_", x)
  row.names(HumanHypoTimepoint) <- row.names(KaZhouAll_ByTimepoint@assays$RNA) 
  HumanHypoAllCells <- merge(HumanHypoAllCells, HumanHypoTimepoint, by = 0)
  row.names(HumanHypoAllCells) <- HumanHypoAllCells$Row.names
  HumanHypoAllCells$Row.names <- NULL
}

BrainSpan <- read.csv("/Users/hannahglover/Downloads/BrainSpan/expression_matrix.csv", head = F, row.names = 1)
BrainSpan_rows <-read.csv("/Users/hannahglover/Downloads/BrainSpan/rows_metadata.csv")
BrainSpan_columns <-read.csv("/Users/hannahglover/Downloads/BrainSpan/columns_metadata.csv")
BrainSpan_columns$SampleID <- paste(BrainSpan_columns$structure_acronym, BrainSpan_columns$age, sep="_")
colnames(BrainSpan) <- BrainSpan_columns$SampleID
BrainSpan$Genes <- BrainSpan_rows$gene_symbol

HumanHypoAllCells$HumanHypo_AllCells_AllTimepoints <- NULL
CombineDatasets <- merge(HumanHypoAllCells, BrainSpan, by.x = 0, by.y = "Genes", all=F)
#CombineDatasets$Dups <- duplicated(CombineDatasets$Row.names) | duplicated(CombineDatasets$Row.names, fromLast=T)
CombineDatasets<- CombineDatasets %>% group_by(Row.names) %>% summarise_all(sum) %>% as.data.frame()
row.names(CombineDatasets) = CombineDatasets$Row.names
CombineDatasets$Row.names <- NULL
write.csv(CombineDatasets, "DevelopingBrainRNAseq_Hypo_BrainSpan_ForPlot1.csv")


#### HUMAN HYPOTHALAMIC NEURONS #### FOR PLOT 2

EdKaZhouHypoNeurons <- readRDS("/Users/hannahglover/Dropbox/Columbia/Brian Hypothalamus/ScienceAdvances_2023_Revision2/CLEAN_SEURAT/EdKaZhouHypoNeurons.rds")

FetalExtrapolated <- read.csv("/Users/hannahglover/Library/CloudStorage/Box-Box/HG2553 Main Folder/Science Advances/FIG2Meta_OCT23.csv", row.names = 1)
FetalExtrapolated <- FetalExtrapolated %>% dplyr::select(AdultFetal_ExtrapolatedNuclei)
EdKaZhouHypoNeurons <- AddMetaData(EdKaZhouHypoNeurons, FetalExtrapolated, "FetalExtrapolated")
EdKaZhouHypoNeurons$NucleiAge <-paste(EdKaZhouHypoNeurons$FetalExtrapolated, EdKaZhouHypoNeurons$SampleAdult, sep="_")
Idents(EdKaZhouHypoNeurons) <- "FetalExtrapolated"
IdentifiedPops <- subset(EdKaZhouHypoNeurons, idents = c("Unassigned", "Fetal"), invert=T)
CheckPopulations <- as.data.frame(table(IdentifiedPops$NucleiAge ))
CheckPopulations20 <- subset(CheckPopulations, Freq >= 20)
Idents(IdentifiedPops) <- "NucleiAge"
IdentifiedPopsClean <- subset(IdentifiedPops, ident = as.character(CheckPopulations20$Var1))
HumanHypoNeurons <- as.data.frame(rowMeans(IdentifiedPopsClean@assays$RNA@counts))
colnames(HumanHypoNeurons) <- "HumanHypo_Neurons_AllTimepoints"
row.names(HumanHypoNeurons) <- row.names(IdentifiedPopsClean@assays$RNA) #36116 genes

for(x in unique(IdentifiedPopsClean@meta.data$NucleiAge)){
  Idents(IdentifiedPopsClean) <- "NucleiAge"
  EdKaZhouHypoNeurons_ByTimepoint <- subset(IdentifiedPopsClean, idents = x)
  
  HumanHypoTimepoint <- as.data.frame(rowMeans(EdKaZhouHypoNeurons_ByTimepoint@assays$RNA@counts))
  colnames(HumanHypoTimepoint) <- paste0("HYPO_", x)
  row.names(HumanHypoTimepoint) <- row.names(EdKaZhouHypoNeurons_ByTimepoint@assays$RNA) 
  HumanHypoNeurons <- merge(HumanHypoNeurons, HumanHypoTimepoint, by = 0)
  row.names(HumanHypoNeurons) <- HumanHypoNeurons$Row.names
  HumanHypoNeurons$Row.names <- NULL
}

HumanHypoNeurons$HumanHypo_Neurons_AllTimepoints <- NULL
write.csv(HumanHypoNeurons, "DevelopingBrainRNAseq_Hypo_ByNucleiAge_ForPlot2.csv")


  
  
#### DOEGE LAB - IN VITRO NEURONS #### PLOT 3
BB_002 = Read10X(data.dir = "~/Library/CloudStorage/Box-Box/HG2553 Main Folder/Bardet-Biedl scRNAseq/LW002_cellranger_analysis_outs/filtered_feature_bc_matrix")
BB_002 =CreateSeuratObject(counts = BB_002, project = "iPSC Cont")
BB_002[["percent.mt"]] <- PercentageFeatureSet(BB_002, pattern = "^MT-")
BB_002_PassedQC = subset(BB_002, percent.mt < 10 & nFeature_RNA > 1000)

InVitro = as.data.frame(rowMeans(BB_002_PassedQC@assays$RNA@layers$counts))
colnames(InVitro) = "DoegeLab_InVitroARC"
row.names(InVitro) = row.names(BB_002_PassedQC@assays$RNA) #33538 genes
write.csv(InVitro, "DoegeLab_InVitroARCNeurons_ForPlot3.csv")

