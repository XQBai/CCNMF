
#' Normalization for single-cell RNA-seq data
#' @param scRNA_raw_matrix the scRNA-seq matrix with raw counts and WGS genes
#' @param gene_locs the genome reference with gene, chr, start, end
#' @return normalized scRNA matrix
#'
#' @export
#'
Normalize_RNA <- function(scRNA_raw_matrix, gene_locs){

  #Step1
  mat <- Filter_genes_below_mean_exp_cutoff(scRNA_raw_matrix, 0.1)
  #Step2
  mat <- Filter_genes(mat, 3)
  #Step3
  mat_normalize <- normalize_counts_by_seq_depth(mat)
  #Step4
  mat_log <- Logtransform(mat_normalize)
  # Step5
  mat_subtract <- Subtract_ref_mean_exp(mat_log)
  # Step6
  mat_bounds <- Max_threshold_bounds(mat_subtract, 3)
  # Step7: smooth data along chromosome with gene windows
  mat_smooth <- Smooth_by_chromosome(mat_bounds, window_length = 101, gene_locs)
  # Step8: Center cells by median
  mat_center <- Center_cells_across_chromosome(mat_smooth, method = 'median')
  # Step9: Adjustment to subtract the relative normal cells
  mat_subtract_adj <- Subtract_ref_mean_exp(mat_center)
  # Step10: Revert the log transformation
  mat_exp <- Invert_logtransform(mat_subtract_adj) + 1
  return(mat_exp)
}

#' Filter genes that expressed fewer than certain number of cells
#' @param mat the input scRNA matrix
#' @param cutoff the percentage of filtering genes
#' @return matrix after filtering genes
#' @export
Filter_genes_below_mean_exp_cutoff <- function(mat, cutoff){

  average_gene <- rowMeans(mat)
  remain_indices <- which(average_gene > cutoff)
  mat <- mat[remain_indices, ]

  return(mat)
}

#' Filter genes with less express cells
#' @param mat the scRNA-seq matrix
#' @param min_num_cells the minmum number of cells
#' @return matrix after filtering genes
#' @export
Filter_genes <- function(mat, min_num_cells){

  total_genes <- dim(mat)[1]
  genes_exp_cn <- apply(mat, 1, function(x){sum(x>0 & !is.na(x))})
  genes_remain = which(genes_exp_cn > min_num_cells)

  if (length(genes_remain) > 0){
    sprintf("%d genes were remained for downstream analysis", length(genes_remain))
    if (length(genes_remain) == 0){
      stop('All genes were removed.')
    }

  }
  mat <- mat[genes_remain, ]
  return(mat)
}

#' Normalizes count data by total sum scaling
#' @param mat the input scRNA matrix
#' @param normalize_factor the scale factor for normalization with depth
#' @return matrix after normaliztion by scale factor
#' @export
normalize_counts_by_seq_depth <- function(mat, normalize_factor = NA){

  # Sum of reads count per cell
  cs = colSums(mat)

  # make fraction of total counts
  data <- sweep(mat, STATS = cs, MARGIN = 2, FUN='/')

  if (is.na(normalize_factor)){
    normalize_factor = median(cs)
  }
  data <- data * normalize_factor
  return(data)
}

#' Log transformation
#' @param data the input matrix
#' @return the log2 normalized matrix
#' @export
Logtransform <- function(data){
  data <- log2(data + 1)
  return(data)
}

#' Invert Log transformation
#' @param data input matrix
#' @return the matrix after invert log transformation
Invert_logtransform <- function(data){
  data <- 2^data - 1
  return(data)
}

#' Subtract average reference
#' @param data the input matrix
#' @param inv_log if perform invert log transformation or not
#' @return the matrix after subtract the average expression
Subtract_ref_mean_exp <- function(data, inv_log = FALSE){

  average_gene <- apply(data, 1, mean)
  data <- sweep(data, STATS = average_gene, MARGIN = 1, FUN='-')

  return(data)
}

#' Apply max/min threshold bounds
#' @param data input matrix
#' @param threshold Set the max and minmum bounds for matrix
#' @return bounded matrix
Max_threshold_bounds <- function(data, threshold){

  data[data > threshold] <- threshold
  data[data < (-1 * threshold) ] <- -1 * threshold

  return(data)
}

#' Smoothing by chromosome through moving average
#' @param obs_data the input matrix
#' @param window_length the number of genes as window for smoothing
#' @return matrix after smoothing
#' @export
smooth_helper <- function(obs_data, window_length){

  orig_obs_data = obs_data

  nas = is.na(obs_data)

  obs_data = obs_data[!nas]

  obs_length <- length(obs_data)
  end_data <- obs_data

  tail_length = (window_length - 1)/2
  if (obs_length >= window_length){
    nas_tmp = is.na(obs_data)
    vals = obs_data[! nas_tmp]

    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))

    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)

    smoothed = stats::filter(vals, custom_filter, sides=2)

    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]

    obs_data[! nas_tmp] = vals
    end_data <- obs_data
  }
  obs_count <- length(obs_data)

  numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))

  # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
  iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))

  for (tail_end in seq_len(iteration_range)) {
    end_tail = obs_count - tail_end + 1

    d_left = tail_end - 1
    d_right = obs_count - tail_end
    d_right = ifelse(d_right > tail_length, tail_length, d_right)

    r_left = tail_length - d_left
    r_right = tail_length - d_right

    denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)

    left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
    right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]

    numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]

    end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
    end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
  }
  # replace the original values by smoothed values

  orig_obs_data[!nas] = end_data

  return(orig_obs_data)
}

#' make the window smooth for matrix
#' @param data the input matrix
#' @param window_length the number of genes for window smooth
#' @return matrix after smooth normaliztion
#' @export
smooth_window <- function(data, window_length){

  if (window_length < 2){
    print('window length < 2, returning orginal unmodified data')
    return(data)
  }

  # Fix ends that couldn't be smoothed
  data_sm <- apply(data, 2, smooth_helper, window_length=window_length)
  # Set back row and columns names
  row.names(data_sm) <- row.names(data)
  colnames(data_sm) <- colnames(data)

  return(data_sm)
}

#' Smooth the matrix on chromosomes
#' @import dplyr
#' @importFrom rlang .data
#' @param data the input matrix
#' @param window_length the number of genes
#' @param gene_locs the geneom reference
#' @param smooth_end if set the smooth end or not
#' @return matrix after smooth on chromosomes
#' @export
Smooth_by_chromosome <- function(data, window_length, gene_locs, smooth_end = TRUE){

  names(gene_locs) <- c('gene', 'chr', 'start', 'end')
  ## Select the genes in matrix
  gene_locs_remain <- gene_locs %>% dplyr::filter(.data$gene %in% rownames(data))
  ## reorder the genes order in matrix
  data <- data[gene_locs_remain$gene, ]

  chrs = unique(gene_locs_remain$chr)
  for (chr in chrs){
    chr_genes_id = which(gene_locs_remain$chr == chr)
    sprintf(paste0('smooth_by_chromosome: chr: ', chr))
    chr_data = data[chr_genes_id, ]

    # consider the only one gene case
    if (is.null(nrow(chr_data))){
      chr_data <- as.matrix(t(chr_data))
    }

    if(nrow(chr_data) > 1){
      smoothed_chr_data = smooth_window(chr_data, window_length)
      data[chr_genes_id, ] <- smoothed_chr_data
    }
  }
  return(data)
}

#' Centering cells after smoothing
#' @param data the input matrix
#' @param method using 'median' or 'mean' for centering
#' @return matrix after centering among cells across chromosomes
#' @export
Center_cells_across_chromosome <- function(data, method = 'median'){

  if (method == 'median'){
    row_median <- apply(data, 2, function(x){median(x, na.rm = TRUE)})
    data <- t(apply(data, 1, '-', row_median))
  }else if (method == 'mean'){
    row_mean <- apply(data, 2, function(x){mean(x, na.rm = TRUE)})
    data <- t(apply(data, 1, '-', row_mean))
  }
  return(data)
}













