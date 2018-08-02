# Methylation Prediction

The goal of this project is find potential regulator of DNA methylation in cancer patients. As a first step, we examine the LAML data set from [Figueroa et al. 2010](https://www.sciencedirect.com/science/article/pii/S1535610809004206) and aims to
- identify subgroups of methylation sites with similar expression pattern in LAML patients from Figueroa dataset 
- find the gene that are strongly associate with the methylation profile of each subgroup of methylation sites

TCGA LAML dataset is used for validation of the findings. 

## Figueroa data set

- [GSE 14468 Expression profiling by array](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14468)
    - platform: GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
    - coverage: 25626 HpaII amplifiable fragments around promoter regions and imprinted regions of the human genome
    - 524 cases of AML

- [GSE 18700 Genome-wide DNA methylation profiling](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18700)
    - platform: GPL6604	HG17_HELP_Promoter
    - 344 cases of AML, 8 normal control
