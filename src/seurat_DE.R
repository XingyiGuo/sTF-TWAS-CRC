# Load library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(DESeq2)

# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ad <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "AD")
# Subset specific cell type
abs_nl_ad <- subset(nl_ad, Cell_Type == "ABS")
# Perform the DEG analysis
Idents(abs_nl_ad) <- "Sample_Classification"
bulk.abs.de <- FindMarkers(object = abs_nl_ad, ident.1 = "NL", ident.2 = "AD", test.use = "DESeq2")
head(bulk.abs.de, n = 15)
write.csv(bulk.abs.de, "../abs_nl_abs_ad_de_results.csv", row.names = TRUE)

# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ad <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "AD")
# Subset specific cell type
gob_nl_ad <- subset(nl_ad, Cell_Type == "GOB")
# Perform the DEG analysis
Idents(gob_nl_ad) <- "Sample_Classification"
bulk.gob.de <- FindMarkers(object = gob_nl_ad, ident.1 = "NL", ident.2 = "AD", test.use = "DESeq2")
head(bulk.gob.de, n = 15)
write.csv(bulk.gob.de, "../gob_nl_gob_ad_de_results.csv", row.names = TRUE)

# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ad <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "AD")
# Subset specific cell type
stm_nl_ad <- subset(nl_ad, Cell_Type == "STM")
# Perform the DEG analysis
Idents(stm_nl_ad) <- "Sample_Classification"
bulk.stm.de <- FindMarkers(object = stm_nl_ad, ident.1 = "NL", ident.2 = "AD", test.use = "DESeq2")
head(bulk.stm.de, n = 15)
write.csv(bulk.stm.de, "../stm_nl_stm_ad_de_results.csv", row.names = TRUE)


# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ser <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "SER")
# Subset specific cell type
abs_nl_ser <- subset(nl_ser, Cell_Type == "ABS")
# Perform the DEG analysis
Idents(abs_nl_ser) <- "Sample_Classification"
bulk.abs.de <- FindMarkers(object = abs_nl_ser, ident.1 = "NL", ident.2 = "SER", test.use = "DESeq2")
head(bulk.abs.de, n = 15)
write.csv(bulk.abs.de, "../abs_nl_abs_ser_de_results.csv", row.names = TRUE)

# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ad <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "SER")
# Subset specific cell type
gob_nl_ser <- subset(nl_ser, Cell_Type == "GOB")
# Perform the DEG analysis
Idents(gob_nl_ser) <- "Sample_Classification"
bulk.gob.de <- FindMarkers(object = gob_nl_ser, ident.1 = "NL", ident.2 = "SER", test.use = "DESeq2")
head(bulk.gob.de, n = 15)
write.csv(bulk.gob.de, "../gob_nl_gob_ser_de_results.csv", row.names = TRUE)

# Read data
seurat_object <- readRDS("../data/For_SC_DE/seurat_object.rds")
# Subset specific disease stage
nl_ad <- subset(seurat_object, Sample_Classification == "NL" | Sample_Classification == "SER")
# Subset specific cell type
stm_nl_ser <- subset(nl_ser, Cell_Type == "STM")
# Perform the DEG analysis
Idents(stm_nl_ser) <- "Sample_Classification"
bulk.stm.de <- FindMarkers(object = stm_nl_ser, ident.1 = "NL", ident.2 = "SER", test.use = "DESeq2")
head(bulk.stm.de, n = 15)
write.csv(bulk.stm.de, "../stm_nl_stm_ser_de_results.csv", row.names = TRUE)
