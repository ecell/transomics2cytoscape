library(transomics2cytoscape)

start.time <- Sys.time()

networkDataDir <- tempfile(); dir.create(networkDataDir)
networkLayers <- "sixfold_yugi2014.tsv"
file.copy(networkLayers, networkDataDir)
stylexml <- system.file("extdata/usecase1", "yugi2014.xml", package = "transomics2cytoscape")
suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
layer1to2 <- "k2e_sixfold.tsv"
layer3to2 <- "sixfold_allosteric_ec2rea.tsv"
file.copy(layer1to2, networkDataDir)
file.copy(layer3to2, networkDataDir)
suid <- createTransomicEdges(suid, layer1to2)
suid <- createTransomicEdges(suid, layer3to2)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
