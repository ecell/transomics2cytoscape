# tests disabled as the Bioconductor build system does not support Cytoscape installation

# test_that("getLayeredNodes and addTransomicsEdges work", {
#     tryCatch({
#         RCy3::cytoscapePing()
#     }, error = function(e) {
#         stop("can't connect to Cytoscape. Please check that Cytoscape is up and running.")
#     })
#     checkCyApps()
#     kinase2enzyme <- system.file('extdata', 'kinase_enzyme.txt',
#                                  package = 'transomics2cytoscape')
#     pathwayZheights <- c(rno00010=1, rno00010=200, rno04910=400, rno04910=600)
#     pathwayIDs <- names(pathwayZheights)
#     lapply(pathwayIDs, transomics2cytoscape:::getKgml)
#     suIDs <- lapply(pathwayIDs, transomics2cytoscape:::importKgml)
#     nodetables <- mapply(transomics2cytoscape:::getNodeTableWithZheight, suIDs, pathwayZheights)
#     edgetables <- lapply(suIDs, transomics2cytoscape:::getEdgeTable)
#     layeredNodes <- transomics2cytoscape:::getLayeredNodes(nodetables)
#     transomicsEdges <- transomics2cytoscape:::addTransomicsEdges(edgetables,
#                                                                 kinase2enzyme, layeredNodes)
#     
#     expect_equal(nrow(layeredNodes), sum(unlist(lapply(nodetables, nrow))))
#     expect_gt(nrow(transomicsEdges), sum(unlist(lapply(edgetables, nrow))))
#     
#     expect_true("id" %in% names(layeredNodes))
#     expect_true("source" %in% names(transomicsEdges))
#     expect_true("target" %in% names(transomicsEdges))
#     
#     expect_true(all(transomicsEdges[, "source"] %in% layeredNodes[, "id"]))
#     expect_true(all(transomicsEdges[, "target"] %in% layeredNodes[, "id"]))
# })
