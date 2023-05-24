#' @description  Merge small bins to large segement
#' @import dplyr
#' @importFrom rlang .data
#' @import stats
#'
#' @param  scdna_matrix: name of copy number variation file
#' @param  bin_size: number of #kb bins
#' @return scdna_matrix_merge  CNV matrix with 1Mb segments
#' @export
#' @examples
#' scdna_object <- Merge_bins2segments(scdna_matrix, bin_size = 50)
Merge_bins2segments <- function(scdna_matrix, bin_size=50){

  genome_reference <- scdna_matrix %>%
    dplyr::select(c('chr', 'start', 'end', 'bin')) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(bin=floor(row_number()/bin_size)) %>%
    dplyr::mutate(bin_corrected = cumsum(c(0,as.numeric(diff(.data$bin))!=0)))

  #aggregating 20kb bins to 1MB segment
  scdna_matrix_merge <- scdna_matrix %>%
    dplyr::mutate(bin_corrected = genome_reference$bin_corrected) %>%
    dplyr::group_by(.data$chr, .data$bin_corrected, .drop = FALSE) %>%
    dplyr::summarise_all(list(median)) %>%
    ungroup() %>%
    dplyr::filter(.data$chr!='chrX' & .data$chr!='chrY') %>%
    dplyr::select(-c('chr', 'start', 'end', 'bin', 'bin_corrected'))

    return(scdna_matrix_merge)
}

#' @importFrom mixtools normalmixEM
#' @import stats
#' @import grDevices
#' @import ggplot2
#'
#' @param scdna_matrix_merge 1Mb-segments copy number matrix
#' @return the index of non-replicating cells
#'
#' @description Identify the replicating/noisy cells by computing standard deviation for cells on 50kb segemental-bins
#' @export
Iden_replicating_cells <- function(scdna_matrix_merge){

  cells_sd <- apply(scdna_matrix_merge, 2, sd)
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

  return(index_segement)
}


#' @description Merge the bins corresponding to the same gene
#' @import dplyr
#' @importFrom rlang .data
#' @param scDNA_matrix The bins times cells CNV matrix
#' @param gene_locs the list of bins index on genes
#' @return the CNV matrix on genes
#'
#' @export
Convert_scDNA_genes_matrix <- function(scdna_matrix, gene_locs){

  if(names(gene_locs) != c('gene', 'chr', 'start', 'end')){
    names(gene_locs) <- c('gene', 'chr', 'start', 'end')
  }

  if(!grepl('chr', gene_locs$chr)){
    gene_locs$chr <- paste0('chr', gene_locs$chr)
  }

  gene_locs <- gene_locs %>% dplyr::filter(.data$chr %in% paste0('chr', seq(1, 22, 1)))

  genome_reference <- scdna_matrix %>%
    dplyr::select(c('chr', 'start', 'end', 'bin')) %>%
    dplyr::filter(!(.data$chr %in% c('chrX', 'chrY')))

  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(genome_reference, keep.extra.columns = TRUE)
  gene_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_locs, keep.extra.columns = TRUE, ignore.strand = TRUE)

  # Find the overlaps between gene and chr regions based on annotation
  olaps <- IRanges::findOverlaps(gene_gr, cnv_gr)

  df_gene <- data_frame(gene = gene_locs$gene[queryHits(olaps)],
                        cnv_bins = cnv_gr$bin[subjectHits(olaps)]) %>%
    dplyr::mutate(bin_corrected = cumsum(c(1,as.numeric(diff(.data$cnv_bins))!=0)))

  scDNA_matrix_gene <- scdna_matrix[df_gene$cnv_bins, ] %>%
    dplyr::select(-c('chr', 'start', 'end')) %>%
    dplyr::mutate(gene = df_gene$gene,
                  bin_corrected = df_gene$bin_corrected) %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise_all(list(mean), na.rm=T) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c('bin', 'bin_corrected'))

  return(scDNA_matrix_gene)
}


