# tests are disabled in the Bioconductor build system
# because it does not support Cytoscape installation

test_that("getLayeredNodes and addTransomicsEdges work", {
    tryCatch({
        RCy3::cytoscapePing()
    }, error = function(e) {
        stop("can't connect to Cytoscape. Please check that Cytoscape is up and running.")
    })
    checkCyApps()
    sif <- system.file("extdata","galFiltered.sif",package="RCy3")
    file.copy(sif, ".")
    networkLayers <- system.file("extdata", "networkLayers.tsv", package = "transomics2cytoscape")
    transomicEdges <- transomicEdges <- system.file("extdata", "transomicEdges.tsv", package = "transomics2cytoscape")
    layerTable <- utils::read.table(networkLayers)
    networkFilePaths <- layerTable$V2
    layers <- lapply(networkFilePaths, transomics2cytoscape:::importLayer)
    networkZheights <- layerTable$V3
    layerIndices <- layerTable$V1
    nodetables <- mapply(transomics2cytoscape:::getNodeTableWithZheight, layers, networkZheights,
                         layerIndices)
    layeredNodes <- transomics2cytoscape:::getLayeredNodes(nodetables)
    
    expect_equal(nrow(layeredNodes), sum(unlist(lapply(nodetables, nrow))))
    
    edgetables <- lapply(layers, transomics2cytoscape:::getEdgeTable)
    et <- dplyr::bind_rows(edgetables)
    et["source"] = as.character(et$source)
    et["target"] = as.character(et$target)
    
    expect_true("id" %in% names(layeredNodes))
    expect_true("source" %in% colnames(et))
    expect_true("target" %in% colnames(et))
    
    suID <- RCy3::createNetworkFromDataFrames(layeredNodes, et)
    createTransomicsEdges(suID, transomicEdges)
    
    newet <- RCy3::getTableColumns(table="edge")
    expect_gt(nrow(newet), nrow(et))
})
