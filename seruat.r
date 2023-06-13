library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(scales)


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
metadata$sample="MCI/AD"
metadata$seq_folder=as.character(metadata$seq_folder)
for(i in 1:nrow(metadata)){
        for(j in 1:nrow(phenotype)){
                if(substr(metadata$seq_folder[i],12,23)==phenotype$RUN_ID[j]){metadata$sample[i]=phenotype$status[j]
                        }
                }
        }

seruat_obj@meta.data=metadata
#save(seruat_obj,"seruat_obj.RDS")


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
        geom_vline(xintercept = 1.1) +
        geom_vline(xintercept = 1.2)
dev.off()


# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = seruat_obj, 
                         subset= (nUMI >= 3000) & 
                           (nGene >= 1200) & 
                           (log10GenesPerUMI > 1.1 & log10GenesPerUMI < 1.2) & 
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

# PCA plot based on cell cycle feature genes
pdf("cell_cycle_PCA.pdf")
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")dev.off
dev.off()

options(future.globals.maxSize = 4000 * 1024^2)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("ctrl", "stim")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.gene, s.features=s.gene)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Check which assays are stored in objects
split_seurat$ctrl@assays


