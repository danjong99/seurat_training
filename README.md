# Raw Data:

The raw data used in this sesseion are collected from sorted CD45+ immune cells of mouse C57BL/6 colon. \
Platform: 10x genomics (3' scRNAseq v2 kits) \
TotalSeq Antibodies: Type A from Biolegend \
**_These data are never used for publication yet. So, please use only within this session._** \
**_If you want to use the data for your own purpose, please let BK know for discussion._**

# Prerequisites:

Install R >= version 3.4: https://www.r-project.org/ \
Install Rstudio Desktop (Free): https://rstudio.com/products/rstudio/download/

# Install Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))\
   install.packages("BiocManager")\
BiocManager::install(version = "3.11")

# Update core packages from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))\
    install.packages("BiocManager")\
BiocManager::install()

# How to Install R Packages:
1.	By Bioconductor – https://www.bioconductor.org/
2.	By install.packages()
3.	From developmental tools – devtools()
4.	From downloaded files

# Required R packages:
Go to Seurat website: https://satijalab.org/seurat/install.html \
Install the following packages: \
(‘Seurat’, ‘dplyr’,’ggplot2’,’patchwork’)

Ready to go!!
