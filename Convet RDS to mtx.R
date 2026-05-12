# ================================
# GSE183276 RDS -> MTX Convert
# ================================

install.packages(c("Matrix", "data.table"))

library(Matrix)
library(data.table)

# Set working directory
setwd("C:/Users/miaso/OneDrive/Desktop/Biomarker from RAW")

# Check files
print(list.files())

# Input file
rds_file <- "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS"

# Check file exists
if (!file.exists(rds_file)) {
  stop("File not found!")
}

# Read RDS
obj <- readRDS(rds_file)

# Show object class
print(class(obj))

# Extract counts matrix
if (inherits(obj, "dgCMatrix") ||
    inherits(obj, "dgTMatrix") ||
    inherits(obj, "dgRMatrix")) {

    counts <- obj

} else if ("Seurat" %in% class(obj)) {

    counts <- obj@assays$RNA@counts

} else {

    stop("Unsupported object type")
}

# Convert to sparse matrix
counts <- as(counts, "dgCMatrix")

# Matrix info
cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n")

# Export Matrix Market
writeMM(counts, "counts.mtx")

# Export genes
fwrite(
  data.table(gene = rownames(counts)),
  "genes.tsv",
  sep = "\t",
  col.names = FALSE
)

# Export barcodes
fwrite(
  data.table(barcode = colnames(counts)),
  "barcodes.tsv",
  sep = "\t",
  col.names = FALSE
)

cat("=================================\n")
cat("Conversion completed!\n")
cat("Generated files:\n")
cat(" - counts.mtx\n")
cat(" - genes.tsv\n")
cat(" - barcodes.tsv\n")
cat("=================================\n")