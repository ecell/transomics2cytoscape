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
# devtools::install_github("ecell/transomics2cytoscape", build_vignettes = TRUE)
```

and also you need to install [Cytoscape](https://cytoscape.org/).

## Example

1. Run Cytoscape Desktop
2. Run R[Studio].
3. Run the following R code.

```R
library(transomics2cytoscape)
kinase2enzyme <- system.file("extdata", "kinase_enzyme.txt",
                             package = "transomics2cytoscape")
create3Dcyjs("rno00010", "rno00010", "rno04910", 1, 200, 400,
             kinase2enzyme, "transomics3D.cyjs")
getwd()
```

4. Import transomics3D.cyjs in the `getwd()` directory with Cytoscsape GUI.
5. Select Cy3D network renderer.

