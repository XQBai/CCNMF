library(CCNMF)
library(data.table)
library(dplyr)
library(ggplot2)


#setwd('./CCNMF/example')

path <- '../data'

# Input data
scDNA <- fread(file.path(path, 'NCI_N87/scDNA/NCI_N87.tsv'))
scDNA_barcode <- read.table(file.path(path, 'NCI_N87/scDNA/NCI-N87_cell_barcodes.txt'))

# The bed file when run cellranger-cnv to make copy number matrix
genome_reference <- fread(file.path(path, 'scDNA_GRCh38_cellranger.bed'),header=F, data.table=F)
names(genome_reference) <- c('chr', 'start', 'end', 'bin')

# Merge CNV matrix with small bins to 50kb-segemental bins as the columns
scdna_object <- Merge_bins2segments(scDNA, genome_reference, bin_size = 50)

# Filter replicating/nosiy cells by computing standard deviation for cells on 50kb segemental-bins
index_remain_cells <- Iden_replicating_cells(scdna_object[[1]])

# import the GRch38 genes bed file
gene_locs <- read.csv(file.path(path, 'gene_locs.bed'), sep='\t', header = F)
colnames(gene_locs) <- c('gene', 'chr', 'start', 'end')
gene_locs <- gene_locs %>% filter(chr != 'X' & chr != 'Y')
gene_locs$chr <- paste0('chr', gene_locs$chr)

# Select genes on signal 50kb-segemental bins
index_signal_segement <- Iden_signal_segments(scdna_object[[1]][,index_remain_cells])
selected_scDNA_genes <- Find_genes_on_signal_segments(scdna_object[[2]], index_signal_segement, gene_locs)

# Align the bins bed file with genes bed file
genes_bin <- Align_bins_genes(gene_locs, genome_reference)

# Convert copy number matrix on bins to genes
# This step costs much time (30~40 mins)
scDNA_genes_matrix <- Convert_scDNA_genes_matrix(scDNA, genes_bin)

# single-cell copy number matrix after filtering noisy cells and slecting signal genes
colnames(scDNA_genes_matrix) <- scDNA_barcode$V1
rownames(scDNA_genes_matrix) <- unique(genes_bin$gene)
scDNA_genes_matrix <- as.matrix(scDNA_genes_matrix)

CNVmatrix <- scDNA_genes_matrix[intersect(unique(genes_bin$gene),selected_scDNA_genes$gene), index_remain_cells]

# Input scRNA matrix
## No InputRNA function
RNAmatrix <- InputRNA(file.path(path, 'NCI_N87/scRNA/'))

# Find common genes between scRNA matrix and scDNA matrix, then as input to CCNMF
CNVmatrix_input <- CNVmatrix[intersect(rownames(CNVmatrix), rownames(RNAmatrix)), ]
RNAmatrix_input <- RNAmatrix[intersect(rownames(CNVmatrix), rownames(RNAmatrix)), ]
saveRDS(CNVmatrix_input, file='./CNVmatrix_input.rds')
saveRDS(RNAmatrix_input, file='./RNAmatrix_input.rds')


### Run CCNMF.R
ncluster = 2
ResultsCCNMF <- Integrate_CCNMF(CNVmatrix_input, RNAmatrix_input, ncluster, initial_parameters = 'same-order',
                                initial_coupling = 'default')

Plot_integrated_figure(CNVmatrix_input, RNAmatrix_input, ncluster, ResultsCCNMF)
### Plot the subclones heatmap for CNV data
Plot_CNV_subclone_heatmap(scDNA_genes_matrix, ResultsCCNMF[[5]], gene_locs, index_remain_cells)



