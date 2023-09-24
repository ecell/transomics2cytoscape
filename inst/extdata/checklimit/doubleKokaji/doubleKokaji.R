library(transomics2cytoscape)

start.time <- Sys.time()
networkDataDir <- tempfile(); dir.create(networkDataDir)
tfs <- system.file("extdata/usecase2", "TFs.sif", package = "transomics2cytoscape")
networkLayers <- "double_kokaji2020.tsv"
stylexml <- system.file("extdata/usecase2", "Kokaji2020styles.xml", package = "transomics2cytoscape")
suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

