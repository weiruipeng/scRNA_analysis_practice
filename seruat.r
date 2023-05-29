



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
for(i in 1:nrow(metadata)){
        for(j in 1:nrow(phenotype)){
                if(metadata$seq_folder[i]==phenotype$RUN_ID[j]){metadata$sample[i]=phenotype$status[j]
                        }
                }
        }

seruat_obj@meta.data=metadata
save(seruat_obj,"seruat_obj.RDS")

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
geom_vline(xintercept = 500)
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
        geom_vline(xintercept = 0.8)
dev.off()
