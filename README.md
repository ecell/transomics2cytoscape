# transomics2cytoscape

## Introduction

To understand transomics datasets, [Yugi et al.
2014](https://pubmed.ncbi.nlm.nih.gov/25131207)
has proposed a method for network visualization
by integrating multiple pathways in 3D space.

The 3D network visualization was created by
[VANTED](https://pubmed.ncbi.nlm.nih.gov/23140568)
and manual operation. transomics2cytoscape automatically
creates the network data for Cytoscape and
[Cy3D](http://apps.cytoscape.org/apps/cy3d) renderer
using
[Cytoscape Automation](https://pubmed.ncbi.nlm.nih.gov/31477170).

## Installation

```{R}
install.packages("devtools")
Sys.setlocale(category = "LC_ALL", locale = "us")
devtools::install_github("ecell/transomics2cytoscape", build_vignettes = FALSE)
```

and also you need to install [Cytoscape](https://cytoscape.org/).

## Example

1. Run Cytoscape Desktop
2. Run R[Studio].
3. Run the following R code. This will import multiple networks and integrate the networks to a 3D space. The networks are layered at the Z coordinate specified by a tsv file [networkLayers.tsv](./inst/extdata/networkLayers.tsv), and are connected by the edges (that is, transomics interaction) between the layers specified in [transomicEdges.tsv](./inst/extdata/transomicEdges.tsv). The 3D network will be styled as specified in [transomics.xml](./inst/extdata/transomics.xml).

```R
library(transomics2cytoscape)
library(dplyr)
sif <- system.file("extdata","galFiltered.sif",package="RCy3")
file.copy(sif, ".")
networkLayers <- system.file("extdata", "networkLayers.tsv",
    package = "transomics2cytoscape")
transomicEdges <- system.file("extdata", "transomicEdges.tsv",
    package = "transomics2cytoscape")
stylexml <- system.file("extdata", "transomics.xml",
    package = "transomics2cytoscape")
create3Dnetwork(networkLayers, transomicEdges, stylexml)
```

Then, you should have a 3D view with layered networks and transomics interactions between them.

![](man/figures/4layers.jpg)
