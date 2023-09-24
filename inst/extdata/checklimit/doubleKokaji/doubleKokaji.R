networkDataDir <- tempfile(); dir.create(networkDataDir)
tfs <- system.file("extdata/usecase2", "TFs.sif", package = "transomics2cytoscape")
networkLayers <- "doubleKokaji.tsv"
stylexml <- system.file("extdata/usecase2", "Kokaji2020styles.xml", package = "transomics2cytoscape")
suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
