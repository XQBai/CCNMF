library(Seurat)
library(optparse)
library(dplyr)
library(viridis)
library(cowplot)
library(future)
# library(CCNMF)
library(optparse)

options(future.globals.maxSize = 6000 * 1024^2)
plan('multiprocess', workers= 18)

setwd('./scRNA_infercnv/')
# Input scRNA matrix
epi_rna <- readRDS('../preprocess_data/seurat_epi.rds')

## Extract tumor epi cells
RNAmatrix <- epi_rna@assays$RNA@counts

# Cell annotation (use cluster# from original sample)
ResultsCCNMF <- readRDS('../ResultsCCNMF.rds')
label = paste0('C', ResultsCCNMF[[6]])
condition <- epi_rna$condition
condition[which(epi_rna$condition == 'tumor')] <- label
cell_annot <- data.frame(label = condition)
cell_annot$barcode <- colnames(epi_rna)

a3 <- RNAmatrix[, which(cell_annot$label == 'normal')]
a1 <- RNAmatrix[, which(cell_annot$label == 'C1')]
a2 <- RNAmatrix[, which(cell_annot$label == 'C2')]
RNAmatrix <-cbind(a3, cbind(a2, a1))

cell_annot <- cell_annot %>% arrange(desc(label))

write.table(cell_annot[,c(2,1)], 'idents.txt', sep="\t", row.names=F, col.names=F)


# import the GRch38 genes bed file
gene_locs <- readRDS('../preprocess_data/genes_locs.rds')
gene_locs <- gene_locs %>% filter(gene %in% rownames(RNAmatrix))

write.table(gene_locs, file = 'hg38_gene_locs.txt', quote=FALSE, sep="\t",
            row.names = FALSE, col.names = FALSE)

counts_matrix <- RNAmatrix[gene_locs$gene, ]
saveRDS(counts_matrix, file = 'counts_matrix.rds')
