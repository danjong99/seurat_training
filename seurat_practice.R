### Please see https://satijalab.org/seurat/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install()


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(1001)

## Load 10x genomics form of data
## These are CD45+ cells from steady state C57BL/6 colon
## with CITEseq antibodies.

cmp.data <- Read10X(data.dir = "./raw_data/")
class(cmp.data)
cmp.data[["Gene Expression"]]
class(cmp.data$`Gene Expression`)
dim(cmp.data$`Gene Expression`)
ncol(cmp.data$`Gene Expression`)
nrow(cmp.data$`Gene Expression`)
colnames(cmp.data$`Gene Expression`)
rownames(cmp.data$`Gene Expression`)

cmp.gene.mtx <- cmp.data$`Gene Expression`

## Data size
dense.size <- object.size(as.matrix(cmp.gene.mtx))
dense.size
sparse.size <- object.size(cmp.gene.mtx)
sparse.size
dense.size/sparse.size

## Creating Seurat object
cmp <- CreateSeuratObject(counts = cmp.gene.mtx, min.cells = 3, min.features = 200, project = "practice")
dim(cmp)
dim(cmp.gene.mtx)

### Data Filtering - regex expression
rownames(cmp)[grep(pattern = "^mt",rownames(cmp))]

cmp[["percent.mt"]] <- PercentageFeatureSet(object = cmp, pattern = "^mt-")
VlnPlot(object = cmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cmp <- subset(x = cmp, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)

plot1 <- FeatureScatter(object = cmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = cmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 #+ plot_layout(ncol = 1)

## RNA expression Normalization
cmp <- NormalizeData(object = cmp, normalization.method = "LogNormalize", scale.factor = 10000)

## Finding Variable Genes
cmp <- FindVariableFeatures(object = cmp, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = cmp), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = cmp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 #+ plot_layout(ncol = 1)

# Data Scaling
all.genes <- rownames(x = cmp)
cmp <- ScaleData(object = cmp, features = all.genes)
cmp <- RunPCA(object = cmp, features = VariableFeatures(object = cmp))
print(cmp[["pca"]])
DimPlot(object = cmp, reduction = "pca", dims = c(1:2))
DimHeatmap(object = cmp, dims = c(1:5), cells = 500, balanced = TRUE)
ElbowPlot(object = cmp)

## Clustering
cmp <- FindNeighbors(object = cmp, dims = 1:16)
cmp <- FindClusters(object = cmp, resolution = 0.4)
head(x = Idents(object = cmp), 5)

## Dimensionality Reduction
cmp <- RunUMAP(object = cmp, dims = 1:16)
DimPlot(object = cmp, reduction = "umap", label = T,
        label.size = 5, pt.size = 0.4) #+ NoLegend() #+ NoAxes()
cmp <- RunTSNE(object = cmp, dims = 1:16)
DimPlot(object = cmp, reduction = "tsne", pt.size = 1,
        label = TRUE, label.size = 5)

cmp.all.markers = FindAllMarkers(object = cmp, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
c1vsc2 <- FindMarkers(object = cmp, ident.1 = 1, ident.2 = 2)

top10 <- cmp.all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cmp, features = top10$gene)

write.table(top10, file = "top10_DGE.txt", sep = '\t')


## http://www.immgen.org/ to find cell types
## refer to review paper https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z

## cell numuber in clusters
table(cmp$seurat_clusters)
## Downsampling
small.cmp = subset(x = cmp, downsample = 100)
DoHeatmap(object = small.cmp, features = top10$gene)
FeaturePlot(object = cmp, features = c("Cd24a"), reduction = 'tsne')

## Renaming clusters
new_cell_ids = c("Bcells","MP1","MP2","Tcells","MP3","Plasmacell1","Stromal1","ILCs","DCs","Stromal2",
  "Stromal3","Stromal4","Plasmacell2")
names(new_cell_ids) <- levels(cmp)
cmp <- RenameIdents(cmp, new_cell_ids)
DimPlot(cmp, label = T)

##### How to get count matrix / scaled matrix
Assays(cmp)
cmp[["RNA"]]@scale.data
mtx = GetAssayData(cmp, assay = "RNA")
class(cmp)
class(mtx)
gene_row_sum = log2(rowSums(as.matrix(mtx)))

x <- rnorm(mean = mean(gene_row_sum),
           sd = sd(gene_row_sum), n = length(gene_row_sum))
par(mfrow = c(1,2))
hist(x, main = "Normal Distribution")
abline(v = quantile(x, probs = seq(0,1,0.05))[20], col = "red")
hist(x = gene_row_sum, labels = FALSE, freq = TRUE,
     xlab = "Avg Gene Exp (log2)", main = "Raw Data Distribution",
     col = "darkmagenta", xlim = c(-5,20))
abline(v = quantile(gene_row_sum, probs = seq(0,1,0.05))[20],
       col = "red")
dev.off()
#Fn = ecdf(x)

hist(x = gene_row_sum, labels = FALSE, freq = TRUE,
     xlab = "Avg Gene Exp (log2)", main = "Raw Data Distribution",
     col = "darkmagenta", xlim = c(-5,20))
abline(v = 13, col = "red")
high.detected.genes <- rownames(mtx)[ gene_row_sum > 13 ]
length(high.detected.genes)
DoHeatmap(object = small.cmp, features = high.detected.genes)

####################################################################################################
### Add CITESeq data to existing Seurat Object
####################################################################################################
class(cmp.data)
names(cmp.data)

cmp.adt.mtx <- cmp.data$`Antibody Capture`
cellid <- colnames(cmp)
cmp.adt.mtx <- cmp.adt.mtx[,cellid]
class(cmp.adt.mtx)
dim(cmp.adt.mtx)
nonZero <- rowSums(as.matrix(cmp.adt.mtx)) > 0 #boolean vector
cmp.adt.mtx <- cmp.adt.mtx[nonZero,]
nrow(cmp.adt.mtx)
plot(x = log2(rowSums(as.matrix(cmp.adt.mtx))), y = apply(as.matrix(cmp.adt.mtx),1,sd),
     xlab = "Exp(log2)", ylab = "S.D.")
abline(v = 7)
nonSmall <- log2(rowSums(as.matrix(cmp.adt.mtx))) > 7
cmp.adt.mtx <- cmp.adt.mtx[nonSmall,]
nrow(cmp.adt.mtx)
rownames(cmp.adt.mtx)
cmp.adt.mtx = cmp.adt.mtx[setdiff(rownames(cmp.adt.mtx), c("CD45.1","CD45.2")),]
nrow(cmp.adt.mtx)

### Add ADT matrix to established Seurat Object
cmp[["ADT"]] <- CreateAssayObject(counts = cmp.adt.mtx)
Assays(cmp)
cmp <- NormalizeData(cmp, assay = "ADT", normalization.method = "CLR")
cmp <- ScaleData(object = cmp, assay = "ADT")
rownames(cmp[["ADT"]])

### Visualization
DefaultAssay(cmp) = "ADT"
FeaturePlot(object = cmp, features = c("Tim4","CD169"),
            min.cutoff = 'q05', max.cutoff = 'q99')
RidgePlot(cmp, features = c("Tim4","CD169"), ncol = 2)

rownames(cmp)
FeatureScatter(object = cmp, feature1 = "CD169", feature2 = "Tim4") + xlim(0,4) +
  geom_hline(yintercept = 1.5) + geom_vline(xintercept = 1.4) + ylim(0,4)
FeatureScatter(object = cmp, feature1 = "RatIgG2akappa", feature2 = "Tim4") + xlim(0,4) +
  geom_hline(yintercept = 1.5) + geom_vline(xintercept = 1.4) + ylim(0,4)
FeatureScatter(object = cmp, feature1 = "CD169", feature2 = "RatIgG2bkappa") + xlim(0,4) +
  geom_hline(yintercept = 1.5) + geom_vline(xintercept = 1.4) + ylim(0,4)

cmp.adt.markers <- FindAllMarkers(cmp, assay = "ADT", only.pos = TRUE)
DoHeatmap(object = cmp, features = unique(cmp.adt.markers$gene), assay = "ADT",
          angle = 90) + NoLegend()
cmp.small <- subset(cmp, downsample = 200)
cmp.adt.markers <- FindAllMarkers(cmp.small, assay = "ADT", only.pos = TRUE)
DoHeatmap(object = cmp.small, features = unique(cmp.adt.markers$gene), assay = "ADT",
          angle = 90) + NoLegend()

saveRDS(cmp, file = "seurat_practice.rds")

#DefaultAssay(cmp) = "ADT"
#cmp <- RunPCA(cmp, features = rownames(cmp), reduction.name = "pcaadt", reduction.key = "pcaadt_", 
#               verbose = FALSE)
#DimPlot(cmp, reduction = "pcaadt")

### Extracting ADT mtx from Seurat Object
#adt.data <- GetAssayData(cmp, slot = "data")
#adt.dist <- dist(t(adt.data)) # distance from t-transformed mtx

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
#cmp[["rnaClusterID"]] <- Idents(cmp)

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
#cmp[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
#cmp[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
#cmp <- FindClusters(cmp, resolution = 0.3, graph.name = "adt_snn")

# We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
#clustering.table <- table(Idents(cmp), cmp$rnaClusterID)
#clustering.table

#DimPlot(cmp, reduction = "tsne_adt", label = T) + NoLegend()
#DimPlot(cmp, reduction = "tsne_adt", group.by = "rnaClusterID", label = T) + NoLegend()
#tsne_rnaClusters <- DimPlot(cmp, reduction = "tsne_adt", group.by = "rnaClusterID") + NoLegend()
#tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
#tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 4)

#tsne_adtClusters <- DimPlot(cmp, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()
#tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
#tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "ident", size = 4)

#wrap_plots(list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

### Be Careful with your results
#cmp.small <- subset(cmp, downsample = 200)
#cmp.adt.markers <- FindAllMarkers(cmp.small, assay = "ADT", only.pos = TRUE)
#DoHeatmap(object = cmp.small, features = unique(cmp.adt.markers$gene), assay = "ADT",
#          angle = 90) + NoLegend()
