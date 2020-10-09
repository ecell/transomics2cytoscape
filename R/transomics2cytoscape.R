##' Import multiple KEGG pathways and integrate the pathways
##' into Cy3D renderer
##'
##' @title Create 3D network view for transomics visualization.
##' @param networkDataDir Path of a directory to put the network files
##' of the second column of networkLayers TSV.
##' @param networkLayers Path of a TSV file with the 3 columns (layer index,
##' the network file name in networkDataDir, Z-height of the network).
##' @param transomicEdges Path of a TSV file with the 5 columns
##' (layer index of a source node,
##' name or KEGG object ID that the source node should have,
##' layer index of a target node,
##' name or KEGG object ID that the target node should have,
##' interaction type).
##' @param stylexml Path of a XML file for Cytoscape style
##' @return A SUID of the 3D network. 
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \donttest{
##' networkDataDir <- tempfile(); dir.create(networkDataDir)
##' sif <- system.file("extdata","galFiltered.sif",package="RCy3")
##' file.copy(sif, networkDataDir)
##' networkLayers <- system.file("extdata", "networkLayers.tsv",
##'     package = "transomics2cytoscape")
##' transomicEdges <- system.file("extdata", "transomicEdges.tsv",
##'     package = "transomics2cytoscape")
##' stylexml <- system.file("extdata", "transomics.xml",
##'     package = "transomics2cytoscape")
##' create3Dnetwork(networkDataDir, networkLayers, transomicEdges, stylexml)
##' }

create3Dnetwork <- function(networkDataDir, networkLayers, transomicEdges,
                            stylexml) {
    owd <- setwd(networkDataDir)
    on.exit(setwd(owd))
    tryCatch({
        RCy3::cytoscapePing()
    }, error = function(e) {
        stop("can't connect to Cytoscape. \n
            Please check that Cytoscape is up and running.")
    })
    checkCyApps()
    layerTable <- utils::read.table(networkLayers)
    networkFilePaths <- layerTable$V2
    layers <- lapply(networkFilePaths, importLayer)
    networkZheights <- layerTable$V3
    layerIndices <- layerTable$V1
    nodetables <- mapply(getNodeTableWithZheight, layers, networkZheights,
                        layerIndices)
    layeredNodes <- getLayeredNodes(nodetables)
    message("Getting edge info. This function is kinda slow...")
    edgetables <- lapply(layers, getEdgeTable)
    message("Finished getting edge info! Creating 3D network...")
    et <- dplyr::bind_rows(edgetables)
    et["source"] = as.character(et$source)
    et["target"] = as.character(et$target)
    RCy3::commandsPOST('cy3d set renderer')
    suID <- RCy3::createNetworkFromDataFrames(layeredNodes, et)
    setTransomicStyle(stylexml, suID)
    message("Creating transomic edges to the 3D network...")
    createTransomicsEdges(suID, transomicEdges)
    message("Finished create3Dnetwork function!")
}

createTransomicsEdges <- function(suID, transomicEdges) {
    nt = RCy3::getTableColumns(table = "node", network = suID)
    et = RCy3::getTableColumns(table = "edge", network = suID)
    #message(dim(nt))
    transomicTable = utils::read.table(transomicEdges)
    
    for (sourceLayerIndex in unique(transomicTable$V1)) {
        transomicEdges = dplyr::filter(transomicTable, V1 == sourceLayerIndex)
        sourceNodeTable = dplyr::filter(nt,
                                        LAYER_INDEX == sourceLayerIndex)
        for (i in seq_len(nrow(transomicEdges))) {
            transomicEdge = transomicEdges[i,]
            sourceNodeName = as.character(transomicEdge$V2)
            transomicInteraction = transomicEdge$V5
            for (j in seq_len(nrow(sourceNodeTable))) {
                nodeRow = sourceNodeTable[j,]
                if (sourceNodeTable$IS_KEGG[1]) {
                    if (sourceNodeName %in% unlist(nodeRow$KEGG_ID)) {
                        s = nodeRow$SUID
                        targetLayerIndex = transomicEdge$V3
                        targetNodeName = transomicEdge$V4
                        targetNodeTable = dplyr::filter(nt,
                                                LAYER_INDEX==targetLayerIndex)
                        createTransomicsEdge(et, targetNodeTable,
                                        targetNodeName, s, transomicInteraction)
                    }
                } else {
                    if (grepl(sourceNodeName, nodeRow$name)) {
                        s = nodeRow$SUID
                        targetLayerIndex = transomicEdge$V3
                        targetNodeName = transomicEdge$V4
                        targetNodeTable = dplyr::filter(nt,
                                                LAYER_INDEX==targetLayerIndex)
                        createTransomicsEdge(et, targetNodeTable,
                                        targetNodeName, s, transomicInteraction)
                    }
                }
            }
        }
    }
}

createTransomicsEdge <- function(edges, targetNodeTable, targetNodeName,
                                sourceId, transomicInteraction){
    for (k in seq_len(nrow(targetNodeTable))) {
        nodeRow = targetNodeTable[k,]
        if (targetNodeTable$IS_KEGG[1]) {
            if (targetNodeName %in% unlist(nodeRow$KEGG_ID)) {
                targetId = nodeRow$SUID
                RCy3::addCyEdges(c(sourceId, targetId),
                                edgeType = as.character(transomicInteraction))
            }
        } else {
            if (grepl(targetNodeName, nodeRow$name)) {
                targetId = nodeRow$SUID
                RCy3::addCyEdges(c(sourceId, targetId),
                                edgeType = as.character(transomicInteraction))
            }
        }
    }
}

checkCyApps <- function(){
    apps = RCy3::getInstalledApps()
    # checking Apps
    if (length(grep("Cy3D,", apps)) == 0) {
        message("Cy3D is not installed. transomics2cytoscape installs Cy3D.")
        RCy3::installApp("Cy3D")
    }
    if (length(grep("KEGGScape,", apps)) == 0) {
        message("KEGGScape is not installed. ",
                "transomics2cytoscape installs KEGGScape.")
        RCy3::installApp("KEGGscape")
    }
}

importLayer <- function(networkFilePath){
    fileExtension <- tools::file_ext(networkFilePath)
    if (fileExtension %in% c("sif", "gml", "xgmml", "xml")){
        message("Importing ", networkFilePath)
        suID <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                    networkFilePath, sep=""))
        Sys.sleep(3)
        layer <- list("suID" = suID$networks, "isKEGG" = FALSE)
        return(layer)
    } else {
        message("transomics2cytoscape tries to import ", networkFilePath,
                    " as KEGG pathway.")
        getKgml(networkFilePath)
        networkFilePath <- paste(networkFilePath, ".xml", sep = "")
        message("Importing ", networkFilePath)
        suID <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                    networkFilePath, sep=""))
        layer <- list("suID" = suID$networks, "isKEGG" = TRUE)
        Sys.sleep(3)
        RCy3::setVisualStyle("KEGG Style")
        RCy3::fitContent()
        return(layer)
    }
}

getKgml <- function(pathwayID){
    writeLines(KEGGREST::keggGet(pathwayID, option = "kgml"),
                paste(pathwayID, ".xml", sep = ""))
}

getNodeTableWithZheight <- function(suid, zheight, layerIndex){
    nodetable = RCy3::getTableColumns(table = "node", network = suid$suID)
    IS_KEGG = rep(suid$isKEGG, nrow(nodetable))
    KEGG_NODE_Z = as.character(rep(zheight, nrow(nodetable)))
    LAYER_INDEX = rep(layerIndex, nrow(nodetable))
    if (!(suid$isKEGG)) {
        positionTable = RCy3::getNodePosition(network = suid$suID)
        xLocations = as.integer(positionTable[["x_location"]])
        yLocations = as.integer(positionTable[["y_location"]])
        maxX = max(xLocations)
        minX = min(xLocations)
        maxY = max(yLocations)
        minY = min(yLocations)
        rangeX = maxX - minX
        rangeY = maxY - minY
        if (rangeX > rangeY) {
            adjustedXlocations = (xLocations - minX) / maxX * 1000
            adjustedYlocations = (yLocations - minY) / maxX * 1000
        } else {
            adjustedXlocations = (xLocations - minX) / maxY * 1000
            adjustedYlocations = (yLocations - minY) / maxY * 1000
        }
        nodetable$KEGG_NODE_X = as.character(adjustedXlocations)
        nodetable$KEGG_NODE_Y = as.character(adjustedYlocations)
    }
    return(cbind(nodetable, KEGG_NODE_Z, LAYER_INDEX, IS_KEGG))
}

getEdgeTable <- function(suid){
    et = RCy3::getTableColumns(table = "edge", network = suid$suID)
    ei = RCy3::getEdgeInfo(et$SUID, network = suid$suID)
    et["source"] = unlist(lapply(ei, function(x) x$source))
    et["target"] = unlist(lapply(ei, function(x) x$target))
    return(et)
}

getLayeredNodes <- function(nodetables){
    nodetable3d = bind_rows(nodetables)
    layeredNodes <- nodetable3d %>% rename(id = "SUID")
    layeredNodes["id"] = as.character(layeredNodes$id)
    return(layeredNodes)
}

setTransomicStyle <- function(xml, suid){
    RCy3::importVisualStyles(filename = xml)
    RCy3::setVisualStyle("transomics", network = suid)
    message("Set visual style to 'transomics'")
}
