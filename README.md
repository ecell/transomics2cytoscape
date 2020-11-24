# transomics2cytoscape

## Introduction

Visualization of trans-omic networks helps biological interpretation by
illustrating pathways where the signals are transmitted.

To characterize signals that go across multiple omic layers, [Yugi and
colleagues have proposed a method for network visualization](https://pubmed.ncbi.nlm.nih.gov/25131207/)
by stacking multiple 2D pathways in a 3D space.

The 3D network visualization was realized by [VANTED](https://www.cls.uni-konstanz.de/software/vanted/).
However, the visualization relies on time-consuming manual operation.
Here we propose **transomics2cytoscape**, an R package that automatically creates
3D network visualization in combination with
Cytoscape, [Cy3D App](http://apps.cytoscape.org/apps/cy3d), and
[Cytoscape Automation](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1758-4).

## Installation

1. Install R from https://cran.r-project.org/
    1. (For Windows Users only) In addition, please download and install Rtools 4.0 from https://cran.r-project.org/bin/windows/Rtools/ 
2. Install Cytoscape from https://cytoscape.org/
3. Run Cytoscape.
3. Run R and the following commands in the R console.

```{R}
rm(list = ls())
remove.packages("transomics2cytoscape")
install.packages("devtools")
devtools::install_github("ecell/transomics2cytoscape", ref = "devel", upgrade = "always")
transomics2cytoscape::installCyApps()
```

## Example

1. Run Cytoscape (If Cytoscape is already running, you don't need to run it anymore. transomics2cytoscape works only when 1 Cytoscape [window] is up.)
2. Run R.
3. Run the following R code. This will import multiple networks and integrate the networks to a 3D space.

```R
library(transomics2cytoscape)
networkDataDir <- tempfile(); dir.create(networkDataDir)
networkLayers <- system.file("extdata", "yugi2014.tsv", package = "transomics2cytoscape")
stylexml <- system.file("extdata", "transomics.xml", package = "transomics2cytoscape")
suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
```

Next Run the following R code. This will add edges between the network layers.

```
transomicEdges <- system.file("extdata", "allosteric.tsv", package = "transomics2cytoscape")
suid <- createTransomicEdges(suid, transomicEdges)
```

Then, you should have a 3D view with layered networks and transomic
interactions between them.

