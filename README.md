# codabiocom

**codabiocom** is an R package implementing compositional relevance and association methodologies for microbiome studies.  
It provides tools to identify relevant OTUs, compute association indices, and benchmark model performance across simulation settings.

## 📦 Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Cruzalirio/codabiocom")
```

##  Features

- Compute **association index** for compositional microbiome data  
- Identify **top relevant OTUs** under different covariate-adjusted models  
- Simulation utilities to generate microbiome datasets with known structure  
- Tools to assess:
  - True Positive Rate (TPR)
  - False Discovery Rate (FDR)
  - AUC
  - Execution time  
- Methods to compare **relevance stability** across models  
- Helper functions for processing output and generating summary tables

##  Example Workflow

```r
library(codabiocom)

# Load example data
data(HIV)

# Compute association index
Xnp <- model.matrix(y_HIV ~ MSM_HIV)
output1 <- LRRelev(data = x_HIV, sample = rownames(x_HIV), 
      group = y_HIV,taxa = colnames(x_HIV), otus = colnames(x_HIV), 
      cores = 2,X = Xnp, method = "hanley")
      
```


## 📄 Documentation

Full documentation and function reference are available in the `./man/` directory or via:


## Contributing

Contributions are welcome!  
Please open an issue or submit a pull request with improvements.


##  Authors

Nelson Cruz
Ricardo Alberich
Irene García
Arnau Mir
Raquel Fernandez
Universitat de les Illes Balears  
2026
