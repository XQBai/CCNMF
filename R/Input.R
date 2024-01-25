#' Input the gene expression matrix which is from 10X
#' @importFrom Matrix readMM
#' @importFrom utils read.table
#' @importFrom utils read.csv
#' @param pathRNA the path of the standard three files: barcodes.tsv, genes.tsv, matrix.mtx located
#' @param file_gene the name of the gene file. Default is 'genes.tsv'
#' @param file_barcodes the name of the barcode file. Default is 'barcodes.tsv'
#' @param file_matrix the name of the gene counts matrix. Default is 'matrix.mtx'
#'
#' @return The inputed gene expression matrix in which rows are cells, columns are genes.
#' @export
InputRNA <- function(pathRNA, file_gene = 'genes.tsv', file_barcodes = 'barcodes.tsv', file_matrix = 'matrix.mtx'){

  # Parse gene expression data first
  if(is(pathRNA, 'character') & is(file_gene, 'character') & is(file_barcodes, 'character') || is(file_matrix, 'character')){
    RNAmatrix <- readMM(file.path(pathRNA, file_matrix))
    genes <- read.table(file.path(pathRNA, file_gene))
    barcodes <- read.csv(file.path(pathRNA, file_barcodes), header = FALSE)
    RNAmatrix <- as.matrix(RNAmatrix)

    rownames(RNAmatrix) <- genes$V1
    colnames(RNAmatrix) <- barcodes$V1
  }else{
    stop("Input pathRNA must be characters, file_gene, file_barcodes and file_matrix must be 'genes.tsv', 'barcodes.tsv' and 'matrix.mtx' respectively.")
  }
  return(RNAmatrix)
}

#' Input the copy number copy in scDNA-seq data, which is a csv file.
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @param pathDNA the path of the CNV file located
#' @param file_CNV the name of the CNV file. csv file is suggested.
#' @param Verbose the CNV files is cnvcalls of cellrangers: Otherwise, the Verbose is FALSE when the cnv file is .csv file as the format of example data.
#'
#' @return The CNV matrix in which the first four rows are "chr", "start", "end", "width". The rest of the rows are
#' the single cell samples.
#' @export
InputDNA <- function(pathDNA, file_CNV, Verbose = TRUE){
  if(Verbose == FALSE & is(pathDNA, 'character') & is(file_CNV, 'character')){
    CNVmatrix <- read.csv(file.path(pathDNA, file_CNV))
  }else if(Verbose == TRUE & is(pathDNA, 'character') & is(file_CNV, 'character')){
    CNVmatrix <- read.table(file.path(pathDNA, file_CNV))
    List_index <- Handle_string(rownames(CNVmatrix))
    chr <- List_index$chr
    start <- List_index$start
    end <- List_index$end
    width <- end - start
    CNVmatrix$chr <- chr
    CNVmatrix$start <- start
    CNVmatrix$end <- end
    CNVmatrix$width <- width

  }else{
    stop("Input pathDNA must be characters, file_CNV must be the csv file.")
  }
  return(CNVmatrix)
}

#' @title Handle string
#' @description Convert the 10X cnv data to the input of CCNMF
#' @import stringr
#' @param Chr chromosome
#' @export
Handle_string <- function(Chr){

  chr <- matrix(0, nrow = length(Chr), ncol=1)
  start <- matrix(0, nrow = length(Chr), ncol=1)
  end <- matrix(0, nrow = length(Chr), ncol = 1)

  for (i in  1:length(Chr)){
    if (length(str_extract(Chr[i], "\\d+")) == 1){
      chr[i] <- str_extract(Chr[i], "\\d+")
      l <- str_extract_all(str_sub(Chr[i], 8), "\\d+")
      start[i] <- as.integer( l[[1]][1])
      end[i] <- as.integer(l[[1]][2])
    } else if (length(str_extract(Chr[i], "\\d+")) == 2){
      chr[i] <- str_extract(Chr[i], "\\d+")
      l <- str_extract_all(str_sub(Chr[i], 9), "\\d+")
      start[i] <- as.integer( l[[1]][1])
      end[i] <- as.integer(l[[1]][2])
    }else{
      chr[i] <- str_sub(Chr[i], 6, 6)
      l <- str_extract_all(str_sub(Chr[i], 8), "\\d+")
      start[i] <- as.integer( l[[1]][1])
      end[i] <- as.integer(l[[1]][2])
    }
  }
  return(list(chr=chr, start = start, end=end))
}
