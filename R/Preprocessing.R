
#' Input the gene expression matrix which is from 10X
#' @import Matrix
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
    RNAmatrix <- Matrix::readMM(file.path(pathRNA, file_matrix))
    genes <- utils::read.table(file.path(pathRNA, file_gene))
    barcodes <- utils::read.csv(file.path(pathRNA, file_barcodes), header = FALSE)
    RNAmatrix <- as.matrix(RNAmatrix)
    #rownames(RNAmatrix) <- levels(barcodes$V1)
    rownames(RNAmatrix) <- genes$V1
    colnames(RNAmatrix) <- barcodes$V1
  }else{
    stop("Input pathRNA must be characters, file_gene, file_barcodes and file_matrix must be 'genes.tsv', 'barcodes.tsv' and 'matrix.mtx' respectively.")
  }
  return(RNAmatrix)
}

#' Input the copy number copy in scDNA-seq data, which is a csv file.
#' @param pathDNA the path of the CNV file located
#' @param file_CNV the name of the CNV file. csv file is suggested.
#' @param Verbose the CNV files is cnvcalls of cellrangers: Otherwise, the Verbose is FALSE when the cnv file is .csv file as the format of example data.
#'
#' @return The CNV matrix in which the first four rows are "chr", "start", "end", "width". The rest of the rows are
#' the single cell samples.
#' @export
InputDNA <- function(pathDNA, file_CNV, Verbose = TRUE){
  if(Verbose == FALSE & is(pathDNA, 'character') & is(file_CNV, 'character')){
    CNVmatrix <- utils::read.csv(file.path(pathDNA, file_CNV))
  }else if(Verbose == TRUE & is(pathDNA, 'character') & is(file_CNV, 'character')){
    CNVmatrix <- utils::read.table(file.path(pathDNA, file_CNV))
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

#' @description Convert the 10X cnv data to the input of CCNMF
#' @param chromosome
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

#' Quality control the cells in scRNA-seq matrix
#'
#' @param RNAmatrix The gene expression matrix of scRNA-seq
#' @return The number of cells will be less
#' @export
QualityControl <- function(RNAmatrix){
  cells <- apply(RNAmatrix, 2, function(x){
    sum(x >= 1) >= 2000
  })
  RNAmatrix <- RNAmatrix[,cells]
  return(RNAmatrix)
}

#' @description Filter genes which did not express in less than 10% cells
#'
#' @param RNAmatrix The gene expression matrix of scRNA-seq
#' @return The filtered gene expression
#' @export
gene_filter <- function(RNAmatrix){
  Filter <- apply(RNAmatrix, 1, function(x){
    sum(x > 1) > round(dim(RNAmatrix)[2]*0.1)
  })
  gene_filter_name <- rownames(RNAmatrix)[Filter]
  RNAmatrix <- RNAmatrix[gene_filter_name, ]
  return(RNAmatrix)
}


#' Converting region-based copy number to gene based
#' Find the corresponding cnv regions and gene bins according to the reference genome.
#' This function refers as below:
#' https://kieranrcampbell.github.io/clonealign/preparing_copy_number_data.html
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import IRanges
#' @import dplyr
#' @import org.Hs.eg.db
#' @import GenomicRanges
#' @import stringr
#' @import tidyr
#' @import S4Vectors
#' @importFrom rlang .data
#'
#' @param RNAmatrix the gene expression matrix which is a matrix file
#' @param CNVmatrix the copy number varients matrix which is a dataframe file
#' @param reference_name the name of reference annotations, default is 'hg19'
#'
#' @return RNAmatrix_match in which the number of genes equals to the number of cnv regions in
#'  CNVmatrix_match. More details, the genes in RNAmatrix_match is one-to-one correspondence with the
#'  cnv regions in CNVmatrix_match.
#'  @export
Estimate_Coupled_matrix <- function(RNAmatrix, CNVmatrix, reference_name = 'hg19'){

  # utils::globalVariables(c('TxDb.Mmusculus.UCSC.mm9.knownGene', 'TxDb.Mmusculus.UCSC.mm9.knownGene',
  #                          'entrezgene', 'ensembl_gene_id'))

  if(reference_name == 'hg19'){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if(reference_name == 'hg38'){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if(reference_name == 'mm9'){
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
  } else{
    stop("There is no corresponding reference in TxDb... package, the another function was suggested...")
  }

  # load the corresponding gene annotations
  g <- genes(txdb, single.strand.genes.only = FALSE)
  # Convert 1, 2, 3... to chr1, chr2, chr3...
  df_CNV <- dplyr::mutate(CNVmatrix, chr = paste0("chr", chr))
  # Convert g to a GRanges object
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(df_CNV, keep.extra.columns = TRUE)
  # Find overlaps between gene and region based annotation
  olaps <- IRanges::findOverlaps(g, cnv_gr)
  #Convert this into a gene and copy number data frame
  cell_names <- names(mcols(cnv_gr)@listData)
  df_gene <- dplyr::data_frame(entrezgene = names(g)[queryHits(olaps)])
  for (i in 1:length(cell_names)){
    r <- as.data.frame(mcols(cnv_gr)@listData[i])[subjectHits(olaps),]
    r <- as.data.frame(r)
    colnames(r) <- cell_names[i]
    df_gene <- cbind(df_gene, r)
  }

  # Map it on ensembl gene ids
  entrezgene_ensembl_map <- as.list(org.Hs.egENSEMBL)

  entrezgene_ensembl_map <- lapply(entrezgene_ensembl_map, `[`, 1)
  df_gene <- dplyr::filter(df_gene, entrezgene %in% names(entrezgene_ensembl_map)) %>%
    dplyr::mutate(ensembl_gene_id = unlist(entrezgene_ensembl_map[entrezgene])) %>%
    dplyr::select(ensembl_gene_id, colnames(df_gene)) %>%
    drop_na()

  # Retain the only genes that uniquely mapped.
  df_gene <- dplyr::count(df_gene, ensembl_gene_id) %>%
    dplyr::filter(n == 1) %>%
    dplyr::inner_join(df_gene) %>%
      dplyr::select(-n)
      options(warn = -1)

  # Select the ovelap genes with scRNA-seq data
  gene_index <- rownames(RNAmatrix)
  #common_gene <- intersect(gene_index, df_gene$ensembl_gene_id)
  CNVmatrix_match <- dplyr::filter(df_gene, ensembl_gene_id %in% gene_index)
  CNVmatrix_match1 <- as.matrix(CNVmatrix_match[, 3:dim(CNVmatrix_match)[2]])
  rownames(CNVmatrix_match1) <- CNVmatrix_match$ensembl_gene_id
  if(str_sub(rownames(RNAmatrix)[1], 1, 4) == 'ENSG'){
    RNAmatrix_match <- RNAmatrix[CNVmatrix_match$ensembl_gene_id, ]
  }else{
    stop("Convert the gene names to ensembl_gene_id, for example using getBM function in biomaRt package")
  }
    # stop(getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol", values = rownames(RNAmatrix), mart=ensembl))
   #CoupledMatrix <- as(as.matrix(diag(length(CNVmatrix_match$ensembl_gene_id))), 'dgCMatrix')
  CoupledMatrix <- as.matrix(diag(length(CNVmatrix_match$ensembl_gene_id)))
  return(list(CNVmatrix_match1, RNAmatrix_match, CoupledMatrix))
}

#' Preprocessed scRNA matrix by Seurat including Normalization, scaling and selcting high variable genes.
#' @description Input the path where the 10X RNA-seq located
#' @import Seurat
#' @param pathRNA the folder of scRNA-seq
#' @export
process_RNA_matrix <- function(pathRNA){

  RNAdata <- Read10X(data.dir = pathRNA)
  RNAObj <-Seurat::CreateSeuratObject(counts = RNAdata, project = 'RNAObj', min.cells = round(dim(RNAdata)[2]*0.1), min.features = 2000)
  RNAObj <- Seurat::NormalizeData(RNAObj, normalization.method = 'LogNormalize', scale.factor = 10000)
  RNAObj <- Seurat::FindVariableFeatures(RNAObj, selection.method = 'vst', nfeatures = 2000)
  all.genes <- rownames(RNAObj)
  RNAObj <- Seurat::ScaleData(RNAObj, features = all.genes)
  RNAmatrix <- RNAObj@assays$RNA@scale.data[RNAObj@assays$RNA@var.features, ]
  return(RNAmatrix)
}

#'
#' @import data.table
#' @param pathRNA the folder of RNA
#' @param RNAobject the seurat object of RNA
#' @export
AlignRNAgenes <- function(pathRNA, RNAobject){
  Gene <- utils::read.table(file.path(pathRNA, 'genes.tsv'))
  RNAscale <- RNAobject@assays$RNA@scale.data[RNAobject@assays$RNA@var.features, ]
  RNAscale <- RNAscale[intersect(Gene$V1, rownames(RNAscale)), ]
  #rownames(RNAscale) <- ConvertGenenames(rownames(RNAscale), Gene, Logic = FALSE)
  return(RNAscale)
}

#' @description find corresponding Hug genes and Ens genes
#' @param inputgene input genes
#' @param Gene the genes
#' @param Logic If the index is TRUE, means translate hug genes to Ens genes. If the index is FALSE, means translate Ens to Hug.
#' @export
ConvertGenenames <- function(inputgene, Gene, Logic=TRUE){
  index <- matrix(0, 1, length(inputgene))
  if (Logic == TRUE){
    for (i in 1:length(index)){
      index[i] <- which(Gene$V2 == inputgene[i])
      Genename <- Gene$V1[index]}
  }else{
    for (i in 1:length(index)){
      index[i] <- which(Gene$V1 == inputgene[i])
      Genename <- Gene$V2[index]
    }
  }
  return(Genename)
}

#' @description  Preprocessing single-cell RNA-seq data by Seurat including normalization, scale, select high variable features.
#' @import Seurat
#' @param RNAdata gene expression data
#' @param min.cells the number of cells
#' @param min.features the number of genes
#' @export
run_Seurat_RNA <- function(RNAdata, min.cells = 6, min.features = 0){
  #RNAdata <- Read10X(data.dir = path)
  RNAobject <- Seurat::CreateSeuratObject(counts = RNAdata, project = 'RNAobject', min.cells = min.cells, min.features = min.features)
  RNAobject <- Seurat::NormalizeData(RNAobject, normalization.method = 'LogNormalize', scale.factor = 10000)
  RNAobject <- Seurat::FindVariableFeatures(RNAobject, selection.method = 'vst', nfeatures = 2000)
  all.genes <- rownames(RNAobject)
  RNAobject <- Seurat::ScaleData(RNAobject, features = all.genes)

  RNAobject <- Seurat::RunPCA(RNAobject, features = VariableFeatures(object = RNAobject))
  RNAobject <- Seurat::FindNeighbors(RNAobject, dims = 1:10)
  RNAobject <- Seurat::FindClusters(RNAobject, resolution = 0.2)
  #L <- Idents(RNAobject)
  return(RNAobject)
}




