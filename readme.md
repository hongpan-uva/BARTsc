## Installation and initialization

Install necessary dependencies

```R
install.packages("devtools")

devtools::install_github("immunogenomics/presto")

```

Install BARTsc

```R
devtools::install_github("hongpan-uva/BARTsc")
```

After BARTsc is installed, **for the first time** the user imports it, BARTsc needs to be initialized with function `initialize()`. This step will automatically create a python virtual environment and install the BART2 python module and related dependencies. The user can specify the path for storing relevant data library (recommended) and the path for the module. 

```R
library("BARTsc")

initialize()
```

> Start installing bart2 and related data library? (Yes/no/cancel) Yes
> Installation started...
> Specify the path to store data library 13.3GB (skip to store under bartsc R package directory): /path/of/data/library
> ...
> Specify the path to install bart2 (skip to install under bartsc R package directory): /path/of/module

## Vignette

BARTsc is made to identify functional transcriptional regulators and reveal relative TF activity between cell types in scRNA-seq, scATAC-seq and single-cell multiome data. Click following links to learn how to apply BARTsc analysis on your data.

[Vignette for scRNA-seq](vignettes/scRNA-seq.md)

[Vignette for scATAC-seq](vignettes/scATAC-seq.md)

[Vignette for single-cell multiome (GEX + ATAC)](vignettes/scMultiome.md)
