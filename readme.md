# BARTsc: a comprehensive toolkit for transcriptional regulator analysis of single-cell omics data

BARTsc is an R package that analyzes transcriptional regulators (TRs) based on single-cell omics data. It supports the analysis of scRNA-seq, scATAC-seq and scMultiome (GEX+ATAC) data. As a successor to BART [BART](https://github.com/zang-lab/bart2), which focuses on bulk data analysis, BARTsc learns binding events from a large collection of public ChIP-seq data instead of using TF binding motifs. BARTsc is a great tool for the identification of (1) TRs that contribute to cell cluster signatures, (2) relative activity of TRs across different cell clusters, and (3) key regulators for each cell cluster.



## Installation and initialization

Install necessary dependencies

```R
install.packages("remotes")

remotes::install_github("immunogenomics/presto")
```

Install BARTsc

```R
remotes::install_github("hongpan-uva/BARTsc")
```

After BARTsc is installed, **for the first time** the user imports it, BARTsc needs to be initialized with function `initialize()`. This step will automatically create a python virtual environment and install the BART2 python module and related dependencies. The user can specify the path for storing relevant data library (recommended) and the path for the module. 

```R
library("BARTsc")

initialize()
```

> Start installing bart2 and related data library? (Yes/no/cancel) Yes
> Installation started...
>
> Specify the absolute path to existing data library, e.g. Data/Path/bart2_library (skip if the data library is uninstalled):
>
> Specify the path to store data library 13.3GB (skip to store under bartsc R package directory): /path/of/data/library
> ...
> Specify the path to install bart2 (skip to install under bartsc R package directory): /path/of/module



## Trouble shooting

If you encounter a <font color="red">**“connection timed out”**</font> error while installing the data library, it may indicate that your network is blocking access to the Box server used for hosting the files. 

In this case, you may try 

(1) run `initialize()` with parameter `site="zenodo"`, which allows you to switch to our Zendo mirror. (recommend)

(2) connecting through a VPN and retry the installation. 

(3) download the data library from our OneDrive mirror and install it manually (see below).



## Manual installation of data library (optional)

We recommend you installing data library through `initialize()` command. However, in the cases when there are connection errors, you may want to install it manually.

First, download the data library using any of the following public URLs.

### Public URLs

#### Box

hg38: `https://virginia.box.com/shared/static/2kqczz9gixetcr9p4bl650uyrio5zd33.gz`

mm10: `https://virginia.box.com/shared/static/bxdggnhp4bjz2l5h2zjlisnzp0ac7axf.gz`

#### Zenodo

hg38:`https://zenodo.org/records/18854649/files/hg38_library.tar.gz?download=1`

mm10:`https://zenodo.org/records/18854649/files/mm10_library.tar.gz?download=1`

#### OneDrive

hg38:`https://myuva-my.sharepoint.com/:u:/g/personal/hz9fq_virginia_edu/IQB2IqcSn23wSaVIP9PoUS1iAVgA5x4T06AzsBcrQ0wLiDA?e=RPOgde`

mm10:`https://myuva-my.sharepoint.com/:u:/g/personal/hz9fq_virginia_edu/IQCankHDq3WqQIO7zYFP3EiqAe_bjFK14kWsQ8kJIpYOJZg?e=pO5PK1`

### Decompress the data

Then, run the following code to install the data library.

```bash
# create a target directory
mkdir Data/Path/bart2_library
cd Data/Path/bart2_library

# download data library or, 
# transfer the downloaded data to your target directory
cp OneDrive/Path/hg38_library.tar.gz .
cp OneDrive/Path/mm10_library.tar.gz .

# decompress
tar zxf hg38_library.tar.gz
tar zxf mm10_library.tar.gz

rm hg38_library.tar.gz
rm mm10_library.tar.gz
```

After the library data is set up, when you go through the setup wizard using `initialize()`, input `Data/Path/bart2_library` as existing data library.



## Vignettes

BARTsc is made to identify functional transcriptional regulators and reveal relative TF activity between cell types in scRNA-seq, scATAC-seq and single-cell multiome data. Click following links to learn how to apply BARTsc analysis on your data.

[Vignette for scRNA-seq](vignettes/scRNA-seq.md)

[Vignette for scATAC-seq](vignettes/scATAC-seq.md)

[Vignette for single-cell multiome (GEX + ATAC)](vignettes/scMultiome.md)

