# CCNMF

## Overview

**C**oupled-**C**lone **N**onnegative **M**atrix **F**actorization (CCNMF) is a method for joint inference of clonal structure using paired Single-cell DNA-Seq and RNA-Seq data. The framework is based on optimizing an objective function that maximizes clone structure coherence between single-cell gene expression and copy number profiles, in which the two profiles are copuled by the dosage effect. The coupling dosage effect can be estimated a prior either by a linear regression model using publicly aviable paired RNA and DNA bulk sequencing data ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)), or by using an uninformative prior. The workflow in CCNMF is illustrated in the figure below.
![](https://github.com/XQBai/CCNMF/blob/master/image/CCNMF_flowchart.png)

## Installation

```
install.packages('devtools')
devtools::install_github("XQBai/CCNMF")
```

## Analyze gastric cancer NCI-N87 cell line

* Download NCI-N87 single-cell RNA-seq data from Gene Expression Omnibus, [(GEO accession number GSE142750)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4238683).

* Download NCI-N87 single-cell DNA-seq FASTQ files from National Institute of Health’s SRA repository, [accession number PRJNA498809](https://www.ncbi.nlm.nih.gov/sra/SRX4943580[accn]). Then we utilized Cellranger-DNA pipeline to convert raw FASTQ files to copy number variants
matrix with small non-overlaping bins times cells on the reference genome (GRCh38). The processed scDNA-seq data is available at [CCNMF/data/NCI_N87/scDNA](https://github.com/XQBai/CCNMF/tree/master/data/NCI_N87/scDNA).

* The pipeline using CCNMF to analyzed paired single-cell NCI-N87 cell line datasets is aviailable at [CCNMF/example/NCI_N87_pipeline.R](https://github.com/XQBai/CCNMF/tree/master/example/NCI_N87_pipeline.R).

* Gastric cancer NCI-N87 cell line data was referenced from:

Noemi Andor, Billy T Lau, Claudia Catalanotti, Anuja Sathe, Matthew Kubit, Jiamin Chen, Cristina Blaj, Athena Cherry, Charles D Bangs, Susan M Grimes, Carlos J Suarez, Hanlee P Ji, Joint single cell DNA-seq and RNA-seq of gastric cancer cell lines reveals rules of in vitro evolution, NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa016, https://doi.org/10.1093/nargab/lqaa016

## Examples
### An example of paired scDNA and scRNA simulation:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb)

### An example of CCNMF analysis of simulated paired scRNA and scDNA data

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb)
### An example of CCNMF analysis of real paired scRNA and scDNA data from a cell mixture:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb)

## Reference
Xiangqi Bai, Zhana Duren, Lin Wan and Li C. Xia. [Joint Inference of Clonal Structure using Single-cell Genome and Transcriptome Sequencing Data, bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.04.934455v2)

<!--## License-->
<!--[MIT](https://github.com/XQBai/CCNMF/blob/master/LICENSE)-->
## Contact
### xqbai@amss.ac.cn and li.xia@einsteinmed.org.
