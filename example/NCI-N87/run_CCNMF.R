
library(CCNMF)
library(dplyr)

CNVmatrix_input <- readRDS('CNVmatrix.rds')
RNAmatrix_input <- readRDS('RNAmatrix.rds')

### Run CCNMF.R
ncluster = 3
# user-defined
ResultsCCNMF <- Integrate_CCNMF(CNVmatrix_input, RNAmatrix_input, ncluster,
                                initial_parameters = 'user-defined',
                                lambda1 = 0.01,
                                lambda2 = 1000)
saveRDS(ResultsCCNMF, file = 'ResultsCCNMF.rds')


## Plot tsne
PlotMainResult(CNVmatrix_input, RNAmatrix_input, ResultsCCNMF)

# import the GRch38 genes bed file
scDNA_genes_matrix <- readRDS('scDNA_genes_matrix.rds')
gene_locs <- readRDS('genes_locs.rds')
gene_locs <- gene_locs %>% dplyr::filter(gene %in% rownames(scDNA_genes_matrix))
Plot_CNV_subclone_heatmap(scDNA_genes_matrix, ResultsCCNMF[[5]], gene_locs, is.noise = FALSE)



