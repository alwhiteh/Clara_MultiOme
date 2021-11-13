#Use Renv to install required Bioconductor packages:
# renv::install("bioc::Seurat")
# renv::install("bioc::Signac")
# renv::install("bioc::GenominfoDb")
# renv::install("bioc::BiocManager")
# renv::install("bioc::Signac")
# renv::install("bioc::scran")
# renv::install("bioc::JASPAR2020")
# renv::install("bioc::BSgenome.Mmusculus.UCSC.mm10")
# renv::install("bioc::chromVAR")
# renv::install("bioc::TFBSTools")
# renv::install("bioc::TxDb.Dmelanogaster.UCSC.dm6.ensGene")
# renv::install("bioc::BSgenome.Dmelanogaster.UCSC.dm6")
# renv::install("bioc::org.Dm.eg.db")
# renv::install("bioc::hdf5r")
#renv::install("bioc::sitadela")

# Load packages
library(Seurat) #for RNA
library(hdf5r)
library(Signac) #for ATAC
library(dplyr)
library(Matrix)
library(scran)
library(ggplot2)
library(sctransform)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(renv)
library(BSgenome.Dmelanogaster.UCSC.dm6)
fly_genome <- BSgenome.Dmelanogaster.UCSC.dm6
#library(org.Dm.eg.db)
library(chromVAR)
library(Rcpp)
library(sitadela) #we will use this to get the annotation for Dm6
library(usethis)

#Renv init
#renv::init()
renv::snapshot()
#renv::restore()

# Create Github Repo
#usethis::use_github()

# We will follow the Satija Lab vignette on multiomics from 10x: https://satijalab.org/signac/articles/pbmc_multiomic.html
# Load RNA and ATAC data
counts <- Read10X_h5("/Volumes/TUNEZ/Clara/filtered_feature_bc_matrix.h5")
fragpath <- "/Volumes/TUNEZ/Clara/atac_fragments.tsv.gz"



# get gene annotations for dm6
# We will use sitadela to build our own genome annotation object that includes biotype
buildDir <- file.path(tempdir(),"test_anndb")
dir.create(buildDir)

# The location of the custom database
myDb <- file.path(buildDir,"testann.sqlite")

# Make the annotation
fly_annotation <- loadAnnotation(genome="dm6",refdb="ensembl",type="gene",db=myDb)
fly_annotation # Here we can see what the GRanges object looks like
genome(fly_annotation) <- "dm6"
seqlevelsStyle(fly_annotation) <- "UCSC"
fly_annotation$gene_biotype <- fly_annotation$biotype # This is the variable name that Signac requires

#Create object for RNA Data
fly <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")

# Now we can add the ATAC data to the seurat object
fly[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = fly_annotation
)

DefaultAssay(fly) <- "ATAC"

#calculate nucleosome signal and annotate
fly <- NucleosomeSignal(fly)
Annotation(fly) <- fly_annotation

# Calculate TSS Enrichment
fly <- TSSEnrichment(fly)

# Plot these QC Metrics
VlnPlot(
  object = fly,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# filter out low quality cells
fly2 <- subset(
  x = fly,
  subset = nCount_ATAC < 10000 &
    nCount_RNA < 10000 &
    nCount_ATAC > 10 &
    nCount_RNA > 10 &
    nucleosome_signal < 0.9 #&
    #TSS.enrichment > 1
)
fly2

# call peaks using MACS2
peaks <- CallPeaks(fly2, macs2.path = "/Users/alexanderwhitehead/anaconda3/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_dm6, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(fly2),
  features = peaks,
  cells = colnames(fly2)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
fly2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = fly_annotation
)

# Let's start looking at the RNA data first
DefaultAssay(fly2) <- "RNA"
fly2 <- SCTransform(fly2)
fly2 <- RunPCA(fly2)

# Then find features for the ATAC data
DefaultAssay(fly2) <- "peaks"
fly2 <- FindTopFeatures(fly2, min.cutoff = 5)
fly2 <- RunTFIDF(fly2)
fly2 <- RunSVD(fly2)

DefaultAssay(fly2) <- "SCT"

fly2 <- FindMultiModalNeighbors(
  object = fly2,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight","ATAC.weight"),
  verbose = TRUE)

# build a joint UMAP visualization via RNA 
fly2 <- RunUMAP(
  object = fly2,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)


DimPlot(fly2, label = TRUE, repel = TRUE, reduction = "umap") 


# If we change the default assay, the UMAP projection stays the same :)
DefaultAssay(fly2) <- "peaks"
DimPlot(fly2, label = TRUE, repel = TRUE, reduction = "umap") 
DefaultAssay(fly2) <- "RNA"
DimPlot(fly2, label = TRUE, repel = TRUE, reduction = "umap") 


# ALEX TRIES STUFF HERE::::::
# perform visualization and clustering steps
DefaultAssay(fly2) <- "RNA"
fly2 <- NormalizeData(fly2)
fly2 <- FindVariableFeatures(fly2)
fly2 <- ScaleData(fly2)
fly2 <- RunPCA(fly2, verbose = FALSE)
fly2 <- FindNeighbors(fly2, dims = 1:30)
fly2 <- FindClusters(fly2, resolution = 0.3, verbose = FALSE)
fly2 <- RunUMAP(fly2, dims = 1:30)
DimPlot(fly2, label = TRUE)


#Find Markers for each cluster
fly2.markers <- FindAllMarkers(fly2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
#All_rna2.markers %>% group_by(cluster) %>% slice_max(order_by = "avg_logFC", n = 1) -> Top1Markers
fly2.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2Markers




# first compute the GC content for each peak
# Note to alex : this should already be in the "genes" object !!!
fly2 <- RegionStats(fly2[['peaks']], genome = BSgenome.Dmelanogaster.UCSC.dm6)

# link peaks to genes
fly2 <- LinkPeaks(
  object = fly2,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)
