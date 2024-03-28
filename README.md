# DoddLabTracer
- Author: Zhiwei Zhou (zhouzw@stanford.edu)
- Created: 03/27/2024
- Last modified: 03/27/2024

## Introduction
This workflow is designed to extract the altered metabolite features pairs between unlabeled medium and stable isotope tracer labeled medium.

The **R** based workflow is developed by Zhiwei Zhou. Please feel free to reach out Zhiwei (zhouzw@stanford.edu) if you have any questions.

## Installation
This workflow is based on R, which requires installing some dependent packages first. 

```
# intall public packages
if (!require(devtools)){
    install.packages("devtools")
}

if (!require(BiocManager)){
    install.packages("BiocManager")
}

# Required packages
required_pkgs <- c("crayon", "data.table", "dplyr", "tidyr","readr", "stringr", "tibble", "purrr",
"ggplot2", "ggrepel", "knitr"", "pbapply", "Rdisop", "RaMS", "readxl", "magrittr", "rmarkdown", "writexl")
BiocManager::install(required_pkgs)
devtools::install_github("DoddLab/MassToolsMjhelf")
devtools::install_github("ZhuMetLab/SpectraTools")

# Install DoddLabTracer
devtools::install_github("DoddLab/DoddLabTracer")
```

 
## Example
This is a basic example of HyuA mutants:
```
library(tidyverse)
library(DoddLabTracer)

# analysis of HyuA 
find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
                 peak_table_label = 'hyuA_13CUA_48h_area.txt',
                 path = '~/Project/00_Uric_Acid_project/Data/20240319_isotope_tracing_analysis/hyuA/',
                 control_group = c("WT_UA", "WT_13CUA"),
                 case_group = c('hyuA_UA', 'hyuA_13CUA'),
                 polarity = 'positive',
                 mz_tol = 10,
                 rt_tol = 0.05,
                 p_value_cutoff = 0.05,
                 fold_change_cutoff = 10,
                 p_adjust = FALSE,
                 is_recognize_adducts = TRUE)

```




## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
