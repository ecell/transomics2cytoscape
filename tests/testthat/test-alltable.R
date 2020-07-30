# # tests disabled as the Bioconductor build system does not support Cytoscape installation
# 
# test_that("getLayeredNodes and addTransomicsEdges work", {
#     tryCatch({
#         RCy3::cytoscapePing()
#     }, error = function(e) {
#         stop("can't connect to Cytoscape. Please check that Cytoscape is up and running.")
#     })
#     checkCyApps()
#     sif <- system.file("extdata","galFiltered.sif",package="RCy3")
#     file.copy(sif, ".")
#     networkLayers <- system.file("extdata", "networkLayers.tsv", package = "transomics2cytoscape")
#     transomicEdges <- transomicEdges <- system.file("extdata", "transomicEdges.tsv", package = "transomics2cytoscape")
#     layerTable <- utils::read.table(networkLayers)
#     networkFilePaths <- layerTable$V2
#     layers <- lapply(networkFilePaths, transomics2cytoscape:::importLayer)
#     networkZheights <- layerTable$V3
#     layerIndices <- layerTable$V1
#     nodetables <- mapply(transomics2cytoscape:::getNodeTableWithZheight, layers, networkZheights,
#                          layerIndices)
#     edgetables <- lapply(layers, transomics2cytoscape:::getEdgeTable)
#     layeredNodes <- transomics2cytoscape:::getLayeredNodes(nodetables)
#     allEdges <- transomics2cytoscape:::addTransomicsEdges(edgetables,
#                                 transomicEdges, layeredNodes)
#     allEdges <- apply(allEdges, 2, as.character)
# 
#     expect_equal(nrow(layeredNodes), sum(unlist(lapply(nodetables, nrow))))
#     expect_gt(nrow(allEdges), sum(unlist(lapply(edgetables, nrow))))
# 
#     expect_true("id" %in% names(layeredNodes))
#     expect_true("source" %in% colnames(allEdges))
#     expect_true("target" %in% colnames(allEdges))
# 
#     expect_true(all(allEdges[, "source"] %in% layeredNodes[, "id"]))
#     expect_true(all(allEdges[, "target"] %in% layeredNodes[, "id"]))
# })
