# CCNMF

## Overview

**C**oupled-**C**lone **N**onnegative **M**atrix **F**actorization (CCNMF) is a method for joint inference of clonal structure using paired Single-cell DNA-Seq and RNA-Seq data. The framework is based on optimizing an objective function that maximizes clone structure coherence between single-cell gene expression and copy number profiles, in which the two profiles are coupled by the dosage effect. The coupling dosage effect can be estimated prior either by a linear regression model using publicly available paired RNA and DNA bulk sequencing data ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)), or by using an uninformative prior. The workflow in CCNMF is illustrated in the figure below.
![](https://github.com/XQBai/CCNMF/blob/master/image/CCNMF_flowchart.png)

## Installation

```
install.packages('devtools')
devtools::install_github("labxscut/CCNMF")
```

## Analyze gastric cancer NCI-N87 cell line
### Raw data

* Download NCI-N87 single-cell RNA-seq data from Gene Expression Omnibus [(GSE142750)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4238683).

* Download NCI-N87 single-cell DNA-seq FASTQ files from National Institute of Health’s SRA repository [PRJNA498809](https://www.ncbi.nlm.nih.gov/sra/SRX4943580[accn]).  

* The raw data can be accessed from [NCI_N87 raw data](https://github.com/labxscut/CCNMF/releases/tag/raw_data). 

### Preprocess data 
* The preprocessing pipeline is as [example/NCI-N87/NCI_N87_preprocess.R](https://github.com/labxscut/CCNMF/blob/main/example/NCI-N87/NCI_N87_preprocess.R).

* The processed scDNA-seq and scRNA-seq matrics are available at [CCNMF/data/processed_data/NCI_N87](https://github.com/labxscut/CCNMF/tree/main/data/processed_data/NCI_N87).

### Run CCNMF
* The pipeline using CCNMF to analyze NCI-N87 cell line datasets is available at [CCNMF/example/NCI_N87/Run_CCNMF.R](https://github.com/labxscut/CCNMF/blob/main/example/NCI-N87/run_CCNMF.R).

*** 

## Analyze primary gastric cancer P5931

### Preprocess data

* The processed scDNA-seq and scRNA-seq matrics are available at [CCNMF/data/processed_data/P5931](https://github.com/labxscut/CCNMF/tree/main/data/processed_data/P5931).

### Run CCNMF
* The pipeline using CCNMF to analyze P5931 is aviailable at [CCNMF/example/P5931/Run_CCNMF.R](https://github.com/labxscut/CCNMF/blob/main/example/P5931/run_CCNMF.R).

*** 

## Simulation
### Generate simulated paired scDNA and scRNA

* [https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb)

<!-- 
### Run CCNMF analysis on simulated data
* [https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb)
-->

<!-- 
### An example of CCNMF analysis of real paired scRNA and scDNA data from a cell mixture:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb)
-->

## Reference
Xiangqi Bai, Zhana Duren, Lin Wan, and Li C. Xia. [Joint Inference of Clonal Structure using Single-cell Genome and Transcriptome Sequencing Data](https://academic.oup.com/nargab/article/6/1/lqae017/7606949) ***NAR Genomics and Bioinformatics***

<!--## License-->
<!--[MIT](https://github.com/XQBai/CCNMF/blob/master/LICENSE)-->
## Contact
### xiangqi@stanford.edu and lcxia@scut.edu.cn.
