library(optparse)
library(dplyr)
library(viridis)
library(cowplot)
library(future)
library(infercnv)
library(CCNMF)

options(future.globals.maxSize = 6000 * 1024^2)
plan('multiprocess', workers= 18)

counts_matrix <- readRDS('counts_matrix.rds')
# create the infercnv object 
infercnv_obj = CreateInfercnvObject(counts_matrix, 
                                    annotations_file= "idents.txt", 
                                    delim="\t", 
                                    gene_order_file='hg38_gene_locs.txt', 
                                    ref_group_names="normal")

# perform infercnv operations to reveal cnv signal
infercnv_obj_final = infercnv::run(infercnv_obj, 
                                   cutoff=.1, 
                                   out_dir="test", 
                                   cluster_by_groups=T, 
                                   #num_threads=40, 
                                   denoise=T,
                                   HMM=T)


