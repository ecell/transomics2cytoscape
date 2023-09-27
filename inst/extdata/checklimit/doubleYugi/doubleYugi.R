library(transomics2cytoscape)

start.time <- Sys.time()

networkDataDir <- tempfile(); dir.create(networkDataDir)
networkLayers <- "double_yugi2014.tsv"
stylexml <- system.file("extdata/usecase1", "yugi2014.xml", package = "transomics2cytoscape")
suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
layer1to2 <- "k2e_double.tsv"
suid <- createTransomicEdges(suid, layer1to2)
layer3to2 <- "double_allosteric_ec2rea.tsv"
suid <- createTransomicEdges(layer3to2, "allosteric_ec2rea.tsv")

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
