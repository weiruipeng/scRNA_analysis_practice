library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(RCurl)
library(cowplot)

# refer to https://github.com/hbctraining/scRNA-seq/blob/master/lessons

f = read.csv("../cellranger_SRR.txt",header=F)
file = c()
for(i in 1:nrow(f)){file=c(file,f[i,])}
for(i in file){
        seruat_data <- Read10X(data.dir = paste0("/home/weir/beegfs/scRNA/",i,"/outs/filtered_feature_bc_matrix"))
        seruat_obj <- CreateSeuratObject(counts = seruat_data, min.features = 100,project = file)
        assign(file,seruat_obj)
}

seruat_obj$log10GenesPerUMI <- log10(seruat_obj$nFeature_RNA)/log10(seruat_obj$nCount_RNA)
seruat_obj$mitoRatio <- PercentageFeatureSet(object = seruat_obj, pattern = "^MT-")
seruat_obj$mitoRatio <- seruat_obj@meta.data$mitoRatio/100

phenotype=read.csv("../SRR_pheotype.csv")
metadata=seruat_obj@meta.data
metadata$sample="MCI/AD"
#metadata$seq_folder=as.character(metadata$seq_folder)
metadata$orig.ident=as.character(metadata$orig.ident)
for(i in 1:nrow(metadata)){
        for(j in 1:nrow(phenotype)){
              # if(substr(metadata$seq_folder[i],12,23)==phenotype$RUN_ID[j]){metadata$sample[i]=phenotype$status[j]}
                if(substr(metadata$orig.ident[i],12,23)==phenotype$RUN_ID[j]){metadata$sample[i]=phenotype$status[j]}
                }
        }

metadata <- metadata %>%
  dplyr::rename(
                #seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA
                )

metadata$seq_folder = metadata$orig.ident
metadata$orig.ident = metadata$sample
seruat_obj@meta.data=metadata
Idents(seruat_obj) <- metadata$sample
  seruat_obj@meta.data$orig.ident


# histogram for ctrl count vs. ADI/DD count
pdf("Ncells.pdf")
metadata %>% 
ggplot(aes(x=sample, fill=sample)) + 
geom_bar() +
#heme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold")) +
ggtitle("NCells")
dev.off()

# density plot for the number UMIs/transcripts per cell
pdf("TranscriptPerCell.pdf")
metadata %>% 
ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("log10 Cell density") +
geom_vline(xintercept = 3000)
dev.off()


# Visualize the distribution of genes detected per cell
pdf("GenePerCellhist.pdf")
metadata %>% 
        ggplot(aes(color=sample, x=nGene, fill= sample)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        geom_vline(xintercept = 1200)
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("nUMInGeneRegression.pdf")
metadata %>% 
        ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
        geom_point() + 
        scale_colour_gradient(low = "gray90", high = "black") +
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        geom_vline(xintercept = 1000) +
        geom_hline(yintercept = 500) +
        facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("MitoExpression.pdf")
metadata %>% 
        ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        geom_vline(xintercept = 0.1)
dev.off()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
pdf("ngenesPerUMI.pdf")
metadata %>%
        ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        geom_vline(xintercept = 1) +
        geom_vline(xintercept = 0.8)
dev.off()


# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = seruat_obj, 
                         subset= (nUMI >= 3000) & 
                           (nGene >= 1200) & 
                           (log10GenesPerUMI > 0.8 & log10GenesPerUMI < 10) & 
                           (mitoRatio < 0.10))

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


seurat_phase <- NormalizeData(filtered_seurat)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can segregate this list into markers of G2/M phase and markers of S phase
s.gene=cc.genes$s.genes
g2m.gene=cc.genes$g2m.genes

seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.gene, 
                                 s.features = s.gene)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
Idents(seurat_phase)=metadata$sample
# perform PCA
seurat_phase <- RunPCA(seurat_phase)

# PCA plot based on cell cycle feature genes
pdf("cell_cycle_PCA.pdf")
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

options(future.globals.maxSize = 4000 * 1024^2)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("ctrl", "MCI/AD")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.gene, s.features=s.gene)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Check which assays are stored in objects
split_seurat$ctrl@assays



# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

Idents(seurat_integrated) = metadata$sample
  

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# Plot PCA
pdf("ctrlvsMCIAD_PCA.pdf")
PCAPlot(seurat_integrated,
        split.by = "sample") 
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP            
pdf("ctrlvsMCIAD_UMAP.pdf")
DimPlot(seurat_integrated)        
dev.off()

pdf("heatmap_for_eachPC.pdf")
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

pdf("ElbowPlotPlot.pdf")
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolutions
head(seurat_integrated@meta.data)

Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
# Plot the UMAP
pdf("clusterUMAP_resolutation08.pdf")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
pdf("clusterUMAP_resolutation04.pdf")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()


# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
head(n_cells)


# UMAP of cells in each cluster by sample
pdf("clusterUMAP_status.pdf")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
dev.off()

pdf("clusterUMAP_cellcycle.pdf")
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
dev.off()


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
pdf("metadataUMAP.pdf")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()


# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:4),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
# Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
pdf("UMAPforeachPC.pdf")
map(paste0("PC_", 1:4), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
dev.off()

# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:4, nfeatures = 5)


# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

pdf("feature.pdf")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

# find different expression gene between clusteres
DefaultAssay(seurat_integrated) <- "RNA"
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

cluster01_markers <- FindMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   ident.2 = 1,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

# FindMarkers will find markers between two different identity groups - you have to specify both identity groups. This is useful for comparing the differences between two specific groups.

# FindAllMarkers will find markers differentially expressed in each identity group by comparing it to all of the others - you don't have to manually define anything. Note that markers may bleed over between closely-related groups - they are not forced to be specific to only one group. This is what most people use (and likely what you want).

# FindConservedMarkers will find markers that are conserved between two groups - this can be useful if you want to find markers that are conserved between a treated and untreated condition for a specific cell type or group of cells. It means they are differentially expressed compared to other groups, but have similar expression between the two groups you're actually comparing.


# Plot interesting marker gene expression of above comparsions. 
pdf("dfexpressmarker_heatmap.pdf")
FeaturePlot(object = seurat_integrated, 
            features = c("LTB", "MAL", "CISH", "PIM1", "IL7R","LEF1"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()

pdf("dfexpressmarker_venn.pdf")
VlnPlot(object = seurat_integrated, 
        features = c("LTB", "MAL", "CISH", "PIM1", "IL7R","LEF1"))
dev.off()



# Determine differentiating markers 
diffexpresscells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,10,18))                  

# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Reorder columns and sort by padj      
cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]

cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj) 

# View data
View(cd4_tcells)

