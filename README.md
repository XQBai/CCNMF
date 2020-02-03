# CCNMF
<center>
**C**oupled-**C**lone **N**onnegative **M**atrix **F**actorization (CCNMF) is a method for joint inference of clonal structure using paired tumor Single-cell DNA-Seq and RNA-Seq data. The framework is based on optimizing an objective function that simultaneously maximizes clone structure coherence between single-cell gene expression matrix and copy number variants matrix, in which the two matrices are copuled by a dosage effect matrix linking expression to copy number. The Coupled matrix can be estimated priorly either by a linear regression model using public paired RNA and DNA bulk sequencing data ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)), or by using an uninformative prior as an identity matrix. Finally, CCNMF simultaneously clusters the coressponding clones for gene expression (scRNA) data and copy number variants (scDNA) data from the same tumor. 
</center>

## Installation

```
devtools::install_github("XQBai/CCNMF")
```

## Authors

Xiangqi Bai, Lin Wan and Li C. Xia.

## License
[MIT](https://github.com/XQBai/CCNMF/blob/master/LICENSE)

### If there are any questions, please contact xqbai@amss.ac.cn or xiangqi@stanford.edu.
