# tests are disabled in the Bioconductor build system
# because it does not support Cytoscape installation

test_that("getNodeTableWithZheight and getEdgeTable work", {
    tryCatch({
        RCy3::cytoscapePing()
    }, error = function(e) {
        stop("can't connect to Cytoscape. Please check that Cytoscape is up and running.")
    })
    pathwayID = "rno00010"
    zheight = 100
    layerIndex = "layer1"
    transomics2cytoscape:::getKgml(pathwayID)
    suID = transomics2cytoscape:::importLayer(pathwayID)
    nodeTable = transomics2cytoscape:::getNodeTableWithZheight(suID, zheight, layerIndex)
    edgeTable = transomics2cytoscape:::getEdgeTable(suID)

    expect_true("KEGG_NODE_Z" %in% names(nodeTable))
    expect_true(all(nodeTable[, "KEGG_NODE_Z"] == 100))
    expect_vector(nodeTable[, "KEGG_NODE_Z"], ptype = character(), size = dim(nodeTable)[1])

    expect_true("source" %in% names(edgeTable))
    expect_true("target" %in% names(edgeTable))
    expect_true(all(edgeTable[, "source"] %in% nodeTable[, "SUID"]))
    expect_true(all(edgeTable[, "target"] %in% nodeTable[, "SUID"]))
})
