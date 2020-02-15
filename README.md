# CCNMF

## Overview

**C**oupled-**C**lone **N**onnegative **M**atrix **F**actorization (CCNMF) is a method for joint inference of clonal structure using paired Single-cell DNA-Seq and RNA-Seq data. The framework is based on optimizing an objective function that maximizes clone structure coherence between single-cell gene expression and copy number profiles, in which the two profiles are copuled by the dosage effect. The coupling dosage effect can be estimated a prior either by a linear regression model using publicly aviable paired RNA and DNA bulk sequencing data ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)), or by using an uninformative prior. The workflow in CCNMF is illustrated in the figure below.
![](https://github.com/XQBai/CCNMF/blob/master/image/CCNMFflow.png)

## Installation

```
devtools::install_github("XQBai/CCNMF")
```
## Usage
### A tutorial on simulation procedures can be found in this notebook:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb)

### A tutoral on CCNMF analyze simulated paired scRNA and scDNA data

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb)
### A tutorial on CCNMF analyze cell mixture paired scRNA and scDNA data can be found in this notebook:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb)

## Paper
[Joint Inference of Clonal Structure using Single-cell DNA-Seq and RNA-Seq data, bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.04.934455v1)

## Authors

Xiangqi Bai, Lin Wan and Li C. Xia.

## License
[MIT](https://github.com/XQBai/CCNMF/blob/master/LICENSE)

### If there are any questions, please contact xqbai@amss.ac.cn or xiangqi@stanford.edu.
