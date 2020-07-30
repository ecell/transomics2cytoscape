# transomics2cytoscape

## Introduction

To understand transomics datasets, [Yugi et al. 2014](https://pubmed.ncbi.nlm.nih.gov/25131207)
has proposed a method for network visualization by integrating multiple pathways in 3D space.
The 3D network visualization was created by
[VANTED](https://pubmed.ncbi.nlm.nih.gov/23140568)
and manual operation.

transomics2cytoscape automatically creates the 3D network visualization with
[Cytoscape](https://cytoscape.org/), 
[Cy3D](http://apps.cytoscape.org/apps/cy3d) renderer, and
[Cytoscape Automation](https://pubmed.ncbi.nlm.nih.gov/31477170).

## Installation

```{R}
install.packages("devtools")
devtools::install_github("ecell/transomics2cytoscape")
```

and also you need to install [Cytoscape](https://cytoscape.org/).

## Example

1. Run Cytoscape Desktop
2. Run R[Studio].
3. Run the following R code. This will import multiple networks and integrate the networks to a 3D space. The networks are layered at the Z coordinate specified by a tsv file [networkLayers.tsv](./inst/extdata/networkLayers.tsv), and are connected by the edges (that is, transomics interaction) between the layers specified in [transomicEdges.tsv](./inst/extdata/transomicEdges.tsv). The 3D network will be styled as specified in [transomics.xml](./inst/extdata/transomics.xml).

```R
library(transomics2cytoscape)
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

The format of `networkLayers.tsv`.
The first column is the network layer ID, the second column is the KEGG pathway ID or an arbitrary network file, and the last column is the Z height of the network.
```
layer1	rno04910	600
layer2	galFiltered.sif	400
layer3	rno00010	200
layer4	rno00010	1
```

The format of `transomicEdges.tsv`.
The first and second column is the information about source node of the transomic interaction. The third and fourth column is about the target node.
The second and fourth column is the id of KEGG object or the name of a node in an arbitrary network file.
The last column is the type of the transomic interaction.
```
layer1	rno:84006	layer2	YMR300C	transomicsType1
layer2	YMR300C	layer3	rno:100364062	transomicsType2
layer3	rno:100364062	layer4	rno:100364062	transomicsType3
```

Then, you should have a 3D view with layered networks and transomics interactions between them.
(Note that you need to perform operations such as zooming out or adjusting the camera angle.)

![4layers](man/figures/4layers.jpg)
