# Part 1: Load and Prepare Multiome Data
source("libraries.r")

gc()

# Get the current working directory
current_directory <- getwd()
# Print the current working directory
print(current_directory)
# Set the working directory to the current directory
setwd(current_directory)

gc()

# Load and setup the 10x multiome object
inputdata.10x <- Read10X_h5("./pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
gc()
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
gc()
# Create Seurat object
obj.multi <- CreateSeuratObject(counts = rna_counts)
gc()
# Get % of mitochondrial genes
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")
gc()
# add the ATAC-seq assay
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
gc()
# Get gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
gc()

# Change style to UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
gc()
# File with ATAC per fragment information file
frag.file <- "./pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
gc()
# Add in ATAC-seq data as ChromatinAssay object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
gc()
# Add the ATAC assay to the multiome object
obj.multi[["ATAC"]] <- chrom_assay
gc()
# Filter ATAC data based on QC metrics
obj.multi <- subset(
  x = obj.multi,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# Save the current session and workspace for Part 1
save.image("part1_session.RData")
save.image("part1_workspace.RData")

gc()

# Part 2: Load and Prepare 10x scATAC-seq Query Data
# Load saved sessions from Part 1
load("part1_session.RData")
load("part1_workspace.RData")

# Load and setup the 10x scATAC-seq query
source("libraries.r")

gc()

# Load ATAC dataset
atac_pbmc_data.10x <- Read10X_h5(filename = "./10k_PBMC_ATAC_nextgem_Chromium_X_filtered_peak_bc_matrix.h5")
gc()
fragpath <- "./10k_PBMC_ATAC_nextgem_Chromium_X_fragments.tsv.gz"
gc()
# Get gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
gc()
# Change to UCSC style 
seqlevelsStyle(annotation) <- 'UCSC'
gc()
# Create ChromatinAssay for ATAC data
atac_pbmc_assay <- CreateChromatinAssay(
  counts = atac_pbmc_data.10x,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
gc()
# Requantify query ATAC to have same features as multiome ATAC dataset
requant_multiome_ATAC <- FeatureMatrix(
  fragments = Fragments(atac_pbmc_assay),
  features = granges(obj.multi[['ATAC']]),
  cells = Cells(atac_pbmc_assay)
)
gc()
# Create assay with requantified ATAC data
ATAC_assay <- CreateChromatinAssay(
  counts = requant_multiome_ATAC,
  fragments = fragpath,
  annotation = annotation
)
gc()
# Create Seurat sbject
obj.atac  <- CreateSeuratObject(counts = ATAC_assay,assay = 'ATAC')
obj.atac[['peak.orig']] <- atac_pbmc_assay
obj.atac <- subset(obj.atac, subset = nCount_ATAC < 7e4 & nCount_ATAC > 2000)
gc()

# Save the current session and workspace for Part 2
save.image("part2_session.RData")
save.image("part2_workspace.RData")

gc()

# Part 3: Integrate and Analyze Data
# Load saved sessions from Parts 1 and 2
load("part1_session.RData")
load("part1_workspace.RData")
load("part2_session.RData")
load("part2_workspace.RData")

# Continue with the integration and analysis...

source("libraries.r")

gc()

h5seurat_file <- "./pbmc_multimodal.h5seurat"
# Read the h5seurat file into a Seurat object
obj.rna <- LoadH5Seurat(h5seurat_file)

gc()

# Preprocessing/normalization for all datasets
# normalize multiome RNA
DefaultAssay(obj.multi) <- "RNA"
obj.multi <- SCTransform(obj.multi, verbose = FALSE)
# normalize multiome ATAC
DefaultAssay(obj.multi) <- "ATAC"
obj.multi <- RunTFIDF(obj.multi)
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
# normalize query
obj.atac <- RunTFIDF(obj.atac)

# Map scATAC-seq dataset using bridge integration
dims.atac <- 2:50
dims.rna <- 1:50
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "SCT"
obj.rna.ext <- PrepareBridgeReference(
  reference = obj.rna, bridge = obj.multi,
  reference.reduction = "spca", reference.dims = dims.rna,
  normalization.method = "SCT"
)

bridge.anchor <- FindBridgeTransferAnchors(
  extended.reference = obj.rna.ext, query = obj.atac,
  reduction = "lsiproject", dims = dims.atac)

obj.atac <- MapQuery(
  anchorset = bridge.anchor, reference = obj.rna.ext,
  query = obj.atac,
  refdata = list(
    l1 = "celltype.l1",
    l2 = "celltype.l2",
    l3 = "celltype.l3"),
  reduction.model = "wnn.umap")

dim_plot1 <- DimPlot(
  obj.atac, group.by = "predicted.l2",
  reduction = "ref.umap", label = TRUE
) + ggtitle("ATAC") + NoLegend()

ggsave("dim_plot1.png", plot = dim_plot1, width = 6, height = 4)

# Assessing the mapping
obj.atac <- FindTopFeatures(obj.atac, min.cutoff = "q0")
obj.atac <- RunSVD(obj.atac)
obj.atac <- RunUMAP(obj.atac, reduction = "lsi", dims = 2:50)

dim_plot2 <- DimPlot(obj.atac, group.by = "predicted.l2", reduction = "umap", label = FALSE)
ggsave("dim_plot2.png", plot = dim_plot2, width = 6, height = 4)

coverage_plot1 <- CoveragePlot(
  obj.atac, region  = "PAX5", group.by = "predicted.l1",
  idents = c("B", "CD4 T", "Mono", "NK"), window = 200,
  extend.upstream = -150000)

ggsave("coverage_plot1_PAX5.png", plot = coverage_plot1, width = 6, height = 4)

coverage_plot2 <- CoveragePlot(
  obj.atac, region = "CD8A", group.by = "predicted.l2",
  idents = c("CD8 Naive", "CD4 Naive", "CD4 TCM", "CD8 TCM"),
  extend.downstream = 5000, extend.upstream = 5000)

ggsave("coverage_plot2_CD8A.png", plot = coverage_plot2, width = 6, height = 4)

coverage_plot3 <- CoveragePlot(
  obj.atac, region = "FOXP3", group.by = "predicted.l2",
  idents = c( "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg"),
  extend.downstream = 0, extend.upstream = 0)

ggsave("coverage_plot3_FOXP3.png", plot = coverage_plot3, width = 6, height = 4)

coverage_plot4 <- CoveragePlot(
  obj.atac, region = "RORC", group.by = "predicted.l2",
  idents = c("CD8 Naive", "CD8 TEM", "CD8 TCM", "MAIT"),
  extend.downstream = 5000, extend.upstream = 5000)
  
ggsave("coverage_plot4_RORC.png", plot = coverage_plot4, width = 6, height = 4)

# Save final session and workspace
save.image("final_session.RData")
save.image("final_workspace.RData")
