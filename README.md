# CCNMF

## Overview

**C**oupled-**C**lone **N**onnegative **M**atrix **F**actorization (CCNMF) is a method for joint inference of clonal structure using paired Single-cell DNA-Seq and RNA-Seq data. The framework is based on optimizing an objective function that maximizes clone structure coherence between single-cell gene expression and copy number profiles, in which the two profiles are copuled by the dosage effect. The coupling dosage effect can be estimated a prior either by a linear regression model using publicly aviable paired RNA and DNA bulk sequencing data ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)), or by using an uninformative prior. The workflow in CCNMF is illustrated in the figure below.
![](https://github.com/XQBai/CCNMF/blob/master/image/CCNMFflow.png)

## Installation

```
devtools::install_github("XQBai/CCNMF")
```
## Examples
### An example of paired scDNA and scRNA simulation:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/SimulationProcedures.ipynb)

### An example of CCNMF analysis of simulated paired scRNA and scDNA data

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/CCNMF_analyze_simulated_data.ipynb)
### An example of CCNMF analysis of real paired scRNA and scDNA data from a cell mixture:

[https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb](https://nbviewer.jupyter.org/github/XQBai/CCNMF/blob/master/notebooks/Real_data_analysis.ipynb)

## Reference
[Xiangqi Bai, Lin Wan and Li C. Xia. Joint Inference of Clonal Structure using Single-cell DNA-Seq and RNA-Seq data, bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.04.934455v1)

<!--## License-->
<!--[MIT](https://github.com/XQBai/CCNMF/blob/master/LICENSE)-->
## Contact
### xqbai@amss.ac.cn and li.xia@einsteinmed.org.
