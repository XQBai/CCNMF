
#' @description  Merge small bins to large segement
#' @import dplyr
#' @import stats
#'
#' @param  scdna_matrix: name of copy number variation file
#' @param  genome_reference: name of bed file
#' @param  bin_size: number of #kb bins
#' @return list(scdna_matrix, scdna_matrix_locs) 50kb matrix and gene_locs
#' @export
#' @examples
#' scdna_object <- Merge_bins2segments(scdna_matrix, genome_reference, bin_size = 50)
Merge_bins2segments <- function(scdna_matrix, genome_reference, bin_size=50){

  genome_reference <- genome_reference %>% dplyr::group_by(chr) %>% dplyr::mutate(bin=floor(row_number()/bin_size))
  genome_reference$bin_corrected <- cumsum(c(0,as.numeric(diff(genome_reference$bin))!=0))

  scdna_matrix$chr <- genome_reference$chr
  scdna_matrix$bin <- genome_reference$bin_corrected
  scdna_matrix$start <- genome_reference$start
  scdna_matrix$end <- genome_reference$end
  scdna_matrix$index <- genome_reference$bin

  #coarse graining
  scdna_matrix <- scdna_matrix %>% dplyr::group_by(chr, bin) %>% dplyr::summarise_all(list(median)) %>% dplyr::filter(chr!='chrX' & chr!='chrY')
  scdna_matrix_locs <- genome_reference %>% dplyr::group_by(chr, bin_corrected) %>% dplyr::mutate(start=min(start)) %>% dplyr::mutate(end = max(end))
  scdna_matrix_locs <- scdna_matrix_locs %>%dplyr::group_by(chr, bin_corrected) %>% dplyr::summarise_all(list(median))
  scdna_matrix_locs <- scdna_matrix_locs[, c('chr', 'start', 'end', 'bin')]
  scdna_matrix_locs <- scdna_matrix_locs %>% dplyr::filter(chr != 'chrX' & chr != 'chrY')

  scdna_matrix$chr <- NULL
  scdna_matrix$bin <- NULL
  scdna_matrix$start <- NULL
  scdna_matrix$end <- NULL
  scdna_matrix$index <- NULL
  saveRDS(scdna_matrix, file='scdna_matrix.rds')
  saveRDS(scdna_matrix_locs, file='scdna_matrix_locs.rds')
  return(list(scdna_matrix, scdna_matrix_locs))
}

#' @importFrom mixtools normalmixEM
#' @import stats
#' @import grDevices
#' @import ggplot2
#'
#' @param scdna_matrix 50kb segmental bins copy number matrix
#' @return the index of non-replicating cells
#'
#' @description Identify the replicating/noisy cells by computing standard deviation for cells on 50kb segemental-bins
#' @export
Iden_replicating_cells <- function(scdna_matrix){

  cells_sd <- apply(scdna_matrix, 2, sd)
  ## Model the mixture distribution and removing small group
  estimate_mix_model <- mixtools::normalmixEM(cells_sd)
  index_no_Sphase <- which(estimate_mix_model$lambda == max(estimate_mix_model$lambda))
  index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
  index_remain <- which(index_label == index_no_Sphase)

  pdf('./Cell_standard_deviation_distribution.pdf')
  plot(density(cells_sd), main='Density of sd', ylab = 'Density')
  dev.off()

  pdf('./Cell_standard_deviation_density.pdf')
  mixtools::plot.mixEM(mixtools::normalmixEM(cells_sd), whichplots = 2)
  dev.off()

  return(index_remain)
}

#' @description Identify the signal segments by computing standard deviation for 50kb bins across cells
#' @importFrom mixtools normalmixEM
#' @import stats
#' @import ggplot2
#' @import grDevices
#'
#' @param scdna_matrix the 50kb-segmental bins * non-replicating cells matrix
#' @return The index of selecting signal segments
#'
#' @export
Iden_signal_segments <- function(scdna_matrix){
  scdna_matrix_no_Sphase <- scdna_matrix
  scdna_matrix_no_Sphase[is.na(scdna_matrix_no_Sphase)] <- 0
  arm_sd <- apply(scdna_matrix_no_Sphase, 1, sd)
  pdf('./Segements_sd_distribution.pdf')
  plot(density(arm_sd), main='Density of sd', ylab = 'Density')
  dev.off()

  pdf('./Segements_sd_density.pdf')
  mixtools::plot.mixEM(normalmixEM(arm_sd), whichplots = 2)
  dev.off()

  estimate_mix_model <- normalmixEM(arm_sd)
  index_seg_cluster <- which(estimate_mix_model$lambda == max(estimate_mix_model$lambda))
  index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
  index_segement <- which(index_label != index_seg_cluster)

  index_segement_tmp <- index_segement
  index_seg_cluster_tmp <- index_seg_cluster

  if (length(index_segement) > round(0.2 * dim(scdna_matrix)[1])){
    index_segement <- index_segement
  }else{
    segement_index_remain <- which(index_label == index_seg_cluster_tmp)
    arm_sd_tmp <- arm_sd[segement_index_remain]
    estimate_mix_model_tmp <- normalmixEM(arm_sd_tmp)
    index_seg_cluster_tmp <- which(estimate_mix_model_tmp$lambda == max(estimate_mix_model_tmp$lambda))
    index_ref <- which(estimate_mix_model_tmp$mu == max(estimate_mix_model_tmp$mu))
    if (index_seg_cluster_tmp == index_ref){
      index_segement <- index_segement
    }else{
      index_label_tmp <- apply(estimate_mix_model_tmp$posterior, 1, function(x){which(x == max(x))})
      index_segement_tmp <- segement_index_remain[which(index_label_tmp != index_seg_cluster_tmp)]
      index_segement <- c(index_segement, index_segement_tmp)
      }
  }

  # ## Stop select when length of selected segements is less than 20% of total
  # while(length(index_segement) <= round(0.2 * dim(scdna_matrix)[1])){
  #
  #   if (length(index_segement) > round(0.2 * dim(scdna_matrix)[1])){
  #     break
  #   }
  #   segement_index_remain <- which(index_label == index_seg_cluster_tmp)
  #   arm_sd_tmp <- arm_sd[segement_index_remain]
  #   estimate_mix_model_tmp <- normalmixEM(arm_sd_tmp)
  #   index_seg_cluster_tmp <- which(estimate_mix_model_tmp$lambda == max(estimate_mix_model_tmp$lambda))
  #   index_label_tmp <- apply(estimate_mix_model_tmp$posterior, 1, function(x){which(x == max(x))})
  #   index_segement_tmp <- segement_index_remain[which(index_label_tmp != index_seg_cluster_tmp)]
  #   index_segement <- c(index_segement, index_segement_tmp)
  # }
  return(index_segement)
}

#' @description Select genes on signal 50kb-segemental bins
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @param scdna_matrix_locs the genome reference for segments
#' @param index_segment the index of signal selected segments
#' @param gene_locs the genes reference
#' @return the names of genes on these signal segments
#'
#' @export
Find_genes_on_signal_segments <- function(scdna_matrix_locs, index_segement, gene_locs){

  segements_locs <- scdna_matrix_locs[index_segement, ]
  gene_gr <- GRanges(seqnames =gene_locs$chr, ranges = IRanges(start = gene_locs$start,
                                                                end = gene_locs$end,
                                                                names = gene_locs$gene))
  segements_gr <- GRanges(seqnames =segements_locs$chr, ranges = IRanges(start = segements_locs$start,
                                                                   end = segements_locs$end,
                                                                   names = segements_locs$chr))
  overlaps <- findOverlaps(gene_gr, segements_gr)
  signal_genes <- data.frame(gene =gene_locs$gene[queryHits(overlaps)],
                             chr = segements_locs$chr[subjectHits(overlaps)])
  return(signal_genes)
}

#' @description Align the bins bed file with genes bed file
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rlang .data
#' @import dplyr
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @param gene_locs the genes reference
#' @param genome_reference the genome reference
#' @return maping bins to one gene
#'
#' @export
Align_bins_genes <- function(gene_locs, genome_reference){

  genome_reference <- genome_reference %>% dplyr::filter(.data$chr != 'chrX' & .data$chr != 'chrY')
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(genome_reference, keep.extra.columns = TRUE)
  gene_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_locs, keep.extra.columns = TRUE, ignore.strand = TRUE)

  # Find the overlaps between gene and chr regions based on annotation
  olaps <- IRanges::findOverlaps(gene_gr, cnv_gr)

  df_gene <- data_frame(gene = gene_locs$gene[queryHits(olaps)],
                        cnv_bins = cnv_gr$bin[subjectHits(olaps)])

  unique_genes <- as.matrix(unique(df_gene$gene))
  bin_index <- apply(unique_genes, 1, function(x){df_gene$cnv_bins[which(df_gene$gene == x)]})

  genes_bin <- data_frame(gene=unique_genes, bins = bin_index)
  return(genes_bin)
}

#' @description Merge the bins corresponding to the same gene
#' @import dplyr
#' @importFrom rlang .data
#' @param scDNA The bins times cells CNV matrix
#' @param genes_bin the list of bins index on genes
#' @return the CNV matrix on genes
#'
#' @export
Convert_scDNA_genes_matrix <- function(scDNA, genes_bin){
  scDNA <- as.matrix(scDNA)
  scDNA <- as.data.frame(scDNA[unlist(genes_bin$bins), ])
  tmp_sizes <- lapply(genes_bin$bins,function(x){length(x)})
  bin_index <- rep(1:length(tmp_sizes), unlist(tmp_sizes))
  scDNA$index <- bin_index
  scDNA <- scDNA %>% dplyr::group_by(.data$index) %>% dplyr::summarise_all(list(mean), na.rm=T)
  scDNA$index <- NULL
  saveRDS(scDNA, file='./scDNA_genes_matrix.rds')
  return(scDNA)
}


