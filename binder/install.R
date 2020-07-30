install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("RCy3")
devtools::install_github("ecell/transomics2cytoscape", build_vignettes = FALSE)
