library(data.table)
library(dplyr)
library(ggplot2)
library(CCNMF)

path <- '../data'
# Input data
scDNA <- fread(file.path(path, 'NCI_N87/scDNA/NCI_N87.tsv'))
scDNA_barcode <- read.table(file.path(path, 'NCI_N87/scDNA/NCI-N87_cell_barcodes.txt'))
# The bed file when run cellranger-cnv to make copy number matrix
genome_reference <- fread(file.path(path, 'NCI_N87/scDNA/scDNA_GRCh38_cellranger.bed'), header=F, data.table=F)
names(scDNA) <- scDNA_barcode$V1
scDNA <- scDNA %>% dplyr::mutate(chr = genome_reference$V1,
                                 start = genome_reference$V2, 
                                 end = genome_reference$V3, 
                                 bin = genome_reference$V4)

scdna_matrix_merge <- Merge_bins2segments(scDNA, bin_size = 50)

# Filter replicating/nosiy cells by computing standard deviation for cells on 50kb segemental-bins
index_remain_cells <- Iden_replicating_cells(scdna_matrix_merge)

# import the GRch38 genes bed file
gene_locs <- read.csv(file.path(path, 'gene_locs.bed'), sep='\t', header = F)

scDNA_genes_matrix  <- Convert_scDNA_genes_matrix(scDNA, gene_locs)
gene <- scDNA_genes_matrix$gene
scDNA_genes_matrix$gene <- NULL
scDNA_genes_matrix <- as.matrix(scDNA_genes_matrix)
rownames(scDNA_genes_matrix) <- gene
CNVmatrix <- scDNA_genes_matrix[, index_remain_cells]
saveRDS(CNVmatrix, file = 'scDNA_genes_matrix.rds')

CNVmatrix[which(CNVmatrix > 10)] = 10
selected_genes <- Iden_signal_segments(CNVmatrix)
CNVmatrix <- CNVmatrix[selected_genes, ]

# Input scRNA matrix
path <- '/data/'
raw_count <- InputRNA(file.path(path, 'NCI_N87/scRNA/'))
## Select the cell barcode for downstream analysis
cells_Gphase <- readRDS('/mnt/ix1/Projects/M070_200622_GI_multiomics/Integrated/CCNMF_test/CCNMF_code_2022/NCI-N87/data/RNAsample_G.rds')
raw_count_G <- raw_count[, colnames(cells_Gphase)]
RNA_normalize <- Normalize_RNA(raw_count_G, gene_locs)
selected_genes <- Iden_signal_segments(RNA_normalize)
rna_noise_gene <- setdiff(rownames(RNA_normalize), rownames(RNA_normalize)[selected_genes])
rna_noise_mat <- RNA_normalize[rna_noise_gene, ]

RNA_normalize_up <- (RNA_normalize - mean(rna_noise_mat))/sd(rna_noise_mat)
RNA_normalize_up <- RNA_normalize_up + 2
RNA_normalize_up[which(RNA_normalize_up < 0 )]=0
RNA_normalize_up[which(RNA_normalize_up >= 6)]=6

common_gene <- intersect(rownames(RNA_normalize), rownames(CNVmatrix))
CNVmatrix_input <- CNVmatrix[common_gene, ]
RNAmatrix_input <- as.matrix(RNA_normalize_up[common_gene, ])

saveRDS(CNVmatrix_input, file = 'CNVmatrix.rds')
saveRDS(RNAmatrix_input, file = 'RNAmatrix.rds')

