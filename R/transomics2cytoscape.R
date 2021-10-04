##' Import multiple KEGG pathways and integrate the pathways
##' into Cy3D renderer
##'
##' @title Create 3D network view for transomics visualization.
##' @param networkDataDir Path of a directory to put the network files
##' of the second column of networkLayers TSV.
##' @param networkLayers Path of a TSV file with the 3 columns (layer index,
##' the network file name in networkDataDir, Z-height of the network).
##' @param stylexml Path of a XML file for Cytoscape style
##' @return A SUID of the 3D network. 
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' networkDataDir <- tempfile(); dir.create(networkDataDir)
##' networkLayers <- system.file("extdata", "yugi2014.tsv",
##'     package = "transomics2cytoscape")
##' stylexml <- system.file("extdata", "transomics.xml",
##'     package = "transomics2cytoscape")
##' suid <- create3Dnetwork(networkDataDir, networkLayers, stylexml)
##' }

create3Dnetwork <- function(networkDataDir, networkLayers,
                            stylexml) {
    owd <- setwd(networkDataDir)
    on.exit(setwd(owd))
    tryCatch({
        RCy3::cytoscapePing()
    }, error = function(e) {
        stop("can't connect to Cytoscape. \n
            Please check that Cytoscape is up and running.")
    })
    installCyApps()
    layerTable <- utils::read.table(networkLayers, sep="\t")
    networkSUID = apply(layerTable, 1, importLayer2)
    layerTable <- cbind(layerTable, networkSUID)
    nodetables <- apply(layerTable, 1, getNodeTableWithLayerinfo)
    layeredNodes <- getLayeredNodes(nodetables)
    edgetables <- apply(layerTable, 1, getEdgeTableWithLayerinfo)
    layeredEdges <- getLayeredEdges(edgetables)

    RCy3::commandsPOST('cy3d set renderer')
    suID <- RCy3::createNetworkFromDataFrames(layeredNodes, layeredEdges)
    setTransomicStyle(stylexml, suID)
    return(suID)
}

##' Create Trans-Omic edges between layers of the network
##'
##' @title Create Trans-Omic edges between layers of the network.
##' @param suid A SUID of Cytoscape network
##' @param transomicEdges Path of a TSV file with the 9 columns
##' (layer index of a source node,
##' name or KEGG object ID that the source node should have,
##' layer index of a target node,
##' name or KEGG object ID that the target node should have,
##' interaction type).
##' @return A SUID of the 3D network. 
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' transomicEdges <- system.file("extdata", "allosteric.tsv",
##'     package = "transomics2cytoscape")
##' suid <- createTransomicEdges(suid, transomicEdges)
##' }

createTransomicEdges <- function(suid, transomicEdges) {
    transomicTable <- utils::read.table(transomicEdges, sep="\t")
    nt = RCy3::getTableColumns(table = "node", network = suid)
    et = RCy3::getTableColumns(table = "edge", network = suid)
    addedEdges = list()
    for (i in seq_len(nrow(transomicTable))) {
        row = transomicTable[i,]
        addedEdges = createTransomicEdge(row, nt, et, addedEdges, suid)
    }
    return(suid)
}

createTransomicEdge <- function(row, nt, et, addedEdges, suid) {
    sourceLayerIndex = as.character(row[1])
    sourceTableType = as.character(row[2])
    sourceTableColumnName = as.character(row[3])
    sourceTableValue = as.character(row[4])
    targetLayerIndex = as.character(row[5])
    targetTableType = as.character(row[6])
    targetTableColumnName = as.character(row[7])
    targetTableValue = as.character(row[8])
    transomicEdgeType= as.character(row[9])
    if (sourceTableType == "node" && targetTableType == "edge"){
        addedEdges = createNode2Edge(nt, sourceLayerIndex, sourceTableValue,
                        sourceTableColumnName, et, targetLayerIndex,
                        targetTableValue, targetTableColumnName,
                        transomicEdgeType, addedEdges, suid)
    } else if (sourceTableType == "node" && targetTableType == "node"){
        addedEdges = createNode2Node(nt, sourceLayerIndex, sourceTableValue,
                        sourceTableColumnName, targetLayerIndex,
                        targetTableValue, targetTableColumnName,
                        transomicEdgeType, addedEdges)
    } else if (sourceTableType == "edge" && targetTableType == "edge"){
        addedEdges = createEdge2Edge(et, sourceLayerIndex, sourceTableValue,
                        sourceTableColumnName, targetLayerIndex,
                        targetTableValue, targetTableColumnName,
                        transomicEdgeType, addedEdges, suid)
    }
    return(addedEdges)
}

createNode2Node <- function(nt, sourceLayerIndex, sourceTableValue,
                            sourceTableColumnName, targetLayerIndex,
                            targetTableValue, targetTableColumnName,
                            transomicEdgeType, addedEdges) {
    sourceLayerNt = dplyr::filter(nt, LAYER_INDEX == sourceLayerIndex)
    sourceNodeRows = dplyr::filter(sourceLayerNt, grepl(sourceTableValue,
                                    !!as.name(sourceTableColumnName),
                                    fixed = TRUE))
    targetLayerNt = dplyr::filter(nt, LAYER_INDEX == targetLayerIndex)
    targetNodeRows = dplyr::filter(targetLayerNt, grepl(targetTableValue,
                                    !!as.name(targetTableColumnName),
                                    fixed = TRUE))
    if (nrow(targetNodeRows) > 0) {
        for (i in seq_len(nrow(sourceNodeRows))) {
            sourceSUID = sourceNodeRows[i, 1]
            for (j in seq_len(nrow(targetNodeRows))){
                targetSUID = targetNodeRows[j, 1]
                if (!(list(c(sourceSUID, targetSUID)) %in% addedEdges)) {
                    addedEdges = append(addedEdges,
                                        list(c(sourceSUID, targetSUID)))
                    RCy3::addCyEdges(c(sourceSUID, targetSUID),
                                    edgeType=as.character(transomicEdgeType))
                }
            }
        }
    }
    return(addedEdges)
}

createNode2Edge <- function(nt, sourceLayerIndex, sourceTableValue,
                            sourceTableColumnName, et, targetLayerIndex,
                            targetTableValue, targetTableColumnName,
                            transomicEdgeType, addedEdges, suid){
    layerNt = dplyr::filter(nt, LAYER_INDEX == sourceLayerIndex)
    sourceNodeRows = dplyr::filter(layerNt, grepl(sourceTableValue,
                                    !!as.name(sourceTableColumnName),
                                    fixed = TRUE))
    layerEt = dplyr::filter(et, LAYER_INDEX == targetLayerIndex)
    targetEdgeRows = dplyr::filter(layerEt, grepl(targetTableValue,
                                    !!as.name(targetTableColumnName),
                                    fixed = TRUE))
    if (nrow(targetEdgeRows) > 0) {
        ei = RCy3::getEdgeInfo(unlist(targetEdgeRows["SUID"]))
        midpointNodes = lapply(ei, getMidpointNodeSUID, suid)
        for (i in seq_len(nrow(sourceNodeRows))){
            sourceSUID = sourceNodeRows[i, 1]
            for (j in seq_len(length(midpointNodes))){
                targetSUID = midpointNodes[[j]]
                if (!(list(c(sourceSUID, targetSUID)) %in% addedEdges)) {
                    addedEdges = append(addedEdges,
                                        list(c(sourceSUID, targetSUID)))
                    RCy3::addCyEdges(c(sourceSUID, targetSUID),
                                    edgeType=as.character(transomicEdgeType))
                }
            }
        }
    }
    return(addedEdges)
}

createEdge2Edge <- function(et, sourceLayerIndex, sourceTableValue,
                            sourceTableColumnName, targetLayerIndex,
                            targetTableValue, targetTableColumnName,
                            transomicEdgeType, addedEdges, suid) {
    sourceLayerEt= dplyr::filter(et, LAYER_INDEX == sourceLayerIndex)
    targetLayerEt = dplyr::filter(et, LAYER_INDEX == targetLayerIndex)
    sourceEdgeRows = dplyr::filter(sourceLayerEt, grepl(sourceTableValue,
                                            !!as.name(sourceTableColumnName),
                                            fixed = TRUE))
    targetEdgeRows = dplyr::filter(targetLayerEt, grepl(targetTableValue,
                                            !!as.name(targetTableColumnName),
                                            fixed = TRUE))
    if (nrow(sourceEdgeRows) > 0 && nrow(targetEdgeRows) > 0) {
        sourceEi = RCy3::getEdgeInfo(sourceEdgeRows["SUID"])
        sourceMidpointNodes = lapply(sourceEi, getMidpointNodeSUID, suid)
        targetEi = RCy3::getEdgeInfo(targetEdgeRows["SUID"])
        targetMidpointNodes = lapply(targetEi, getMidpointNodeSUID, suid)
        for (i in seq_len(length(sourceMidpointNodes))){
            sourceSUID = sourceMidpointNodes[[i]]
            for (j in seq_len(length(targetMidpointNodes))) {
                targetSUID = targetMidpointNodes[[j]]
                if (!(list(c(sourceSUID, targetSUID)) %in% addedEdges)) {
                    addedEdges = append(addedEdges,
                                        list(c(sourceSUID, targetSUID)))
                    RCy3::addCyEdges(c(sourceSUID, targetSUID),
                                edgeType=as.character(transomicEdgeType))
                }
            }
        }
    }
    return(addedEdges)
}

getEdgeSourceSUID <- function(edgeInfo){
    sourceNodeSUID = edgeInfo$source
    return(sourceNodeSUID)
}

getMidpointNodeSUID <- function(edgeInfo, suid){
    sourceNodeSUID = edgeInfo$source
    targetNodeSUID = edgeInfo$target
    sourceNodeX = RCy3::cyrestGET(paste('networks',suid,'tables','defaultnode',
                                'rows',sourceNodeSUID,'KEGG_NODE_X',sep="/"))
    targetNodeX = RCy3::cyrestGET(paste('networks',suid,'tables','defaultnode',
                                'rows',targetNodeSUID,'KEGG_NODE_X',sep="/"))
    midNodeX = (as.numeric(sourceNodeX) + as.numeric(targetNodeX)) / 2
    sourceNodeY = RCy3::cyrestGET(paste('networks',suid,'tables','defaultnode',
                                'rows',sourceNodeSUID,'KEGG_NODE_Y',sep="/"))
    targetNodeY = RCy3::cyrestGET(paste('networks',suid,'tables','defaultnode',
                                'rows',targetNodeSUID,'KEGG_NODE_Y',sep="/"))
    midNodeY = (as.numeric(sourceNodeY) + as.numeric(targetNodeY)) / 2
    midNodeZ = RCy3::cyrestGET(paste('networks',suid,'tables','defaultnode',
                                'rows',sourceNodeSUID,'KEGG_NODE_Z',sep="/"))
    n = RCy3::addCyNodes(edgeInfo$SUID)
    midNodeSUID = n[[1]]$SUID
    RCy3::commandsPOST(paste0('node set attribute',
                        ' columnList=KEGG_NODE_X,KEGG_NODE_Y,KEGG_NODE_Z',
                        ' nodeList=SUID:', midNodeSUID,
                        ' valueList=', midNodeX, ',', midNodeY, ',', midNodeZ))
    return(midNodeSUID)
}

##' Convert KEGG enzyme IDs to KEGG reaction IDs
##'
##' @title Convert KEGG enzyme IDs to KEGG reaction IDs.
##' @param tsvFilePath Path of a TSV file with the 9 columns
##' (layer index of a source node,
##' name or KEGG object ID that the source node should have,
##' layer index of a target node,
##' name or KEGG object ID that the target node should have,
##' interaction type).
##' @param columnIndex The column number
##' @param outputFilename The output filename
##' @return None
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' ec = system.file("extdata", "allosteric_ecnumber.tsv",
##'     package = "transomics2cytoscape")
##' ec2reaction(ec, 8, "allosteric.tsv")
##' }

ec2reaction <- function(tsvFilePath, columnIndex, outputFilename) {
    source2target(tsvFilePath, columnIndex, outputFilename, "reaction")
}

ecRow2reaRows <- function(row, columnIndex, ec2rea) {
    ec = row[columnIndex]
    rea = ec2rea[names(ec2rea) == ec]
    rows = do.call("rbind", replicate(length(rea), row, simplify = FALSE))
    return(as.data.frame(cbind(rows, as.vector(rea))))
}

source2target <- function(tsvFilePath, columnIndex, outputFilename, target) {
    transomicTable = utils::read.table(tsvFilePath, sep="\t")
    sourceVec = transomicTable[ , columnIndex]
    sourceVec = unique(sourceVec)
    lastIndex = length(sourceVec)
    if (lastIndex> 100) {
        s2t = c()
        for (i in 1:round(lastIndex/100)) {
            tmp = KEGGREST::keggLink(target, sourceVec[1:100*i])
            s2t = c(s2t, tmp)
        }
        tmp = KEGGREST::keggLink(target,
                            sourceVec[100*round(lastIndex/100)+1:lastIndex])
        s2t = c(s2t, tmp)
    } else {
        s2t = KEGGREST::keggLink(target, sourceVec)
    }
    dflist = apply(transomicTable, 1, srcRow2tgtRows, columnIndex, s2t)
    df = dplyr::bind_rows(dflist)
    df = dplyr::select(df, -one_of("rows"))
    lastColumn = df[ , ncol(df)]
    originalColumn = df[ , columnIndex]
    df[ , columnIndex] = lastColumn
    df[ , ncol(df)] = originalColumn
    utils::write.table(unique(df), file = outputFilename, quote = FALSE,
                sep = '\t', col.names = FALSE, row.names = FALSE)
}

srcRow2tgtRows <- function(row, columnIndex, s2t) {
    srcId = row[columnIndex]
    tgtId = s2t[names(s2t) == srcId]
    rows = do.call("rbind", replicate(length(tgtId), row, simplify = FALSE))
    return(as.data.frame(cbind(rows, as.vector(tgtId))))
}

gene2ec <- function(tsvFilePath, columnIndex, outputFilename) {
    source2target(tsvFilePath, columnIndex, outputFilename, "ec")
}


##' Install the Cytoscape Apps the transomics2cytoscape depends
##'
##' @title Install the Cytoscape Apps the transomics2cytoscape depends.
##' @return None
##' @author Kozo Nishida
##' @export
##' @examples \dontrun{
##' installCyApps()
##' }

installCyApps <- function(){
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

importLayer2 <- function(row){
    fileExtension <- tools::file_ext(row[2])
    if (fileExtension %in% c("sif", "gml", "xgmml", "xml")){
        res <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                        row[2], sep=""))
        Sys.sleep(3)
        nodePositions <- RCy3::getNodePosition(network = res$networks)
        names(nodePositions)[1] <- "KEGG_NODE_X"
        names(nodePositions)[2] <- "KEGG_NODE_Y"
        RCy3::loadTableData(nodePositions, network = res$networks)
        return(res$networks)
    } else {
        getKgml(row[2])
        kgml <- paste(row[2], ".xml", sep = "")
        message("Importing ", kgml)
        res <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                        kgml, sep=""))
        Sys.sleep(3)
        RCy3::setVisualStyle("KEGG Style")
        RCy3::fitContent()
        networkSUID = res$networks
        return(networkSUID)
    }
}

getKgml <- function(pathwayID){
    writeLines(KEGGREST::keggGet(pathwayID, option = "kgml"),
                paste(pathwayID, ".xml", sep = ""))
}

getNodeTableWithLayerinfo <- function(row){
    nt = RCy3::getTableColumns(table = "node", network = as.numeric(row[4]))
    KEGG_NODE_Z = as.character(rep(row[3], nrow(nt)))
    LAYER_INDEX = rep(row[1], nrow(nt))
    return(cbind(nt, KEGG_NODE_Z, LAYER_INDEX))
}

getEdgeTableWithLayerinfo <- function(row){
    et = RCy3::getTableColumns(table = "edge", network = as.numeric(row[4]))
    message("Getting edge info. This function is kinda slow...")
    ei = RCy3::getEdgeInfo(et$SUID, as.numeric(row[4]))
    message("Finished getting edge info.")
    et["source"] = unlist(lapply(ei, function(x) x$source))
    et["target"] = unlist(lapply(ei, function(x) x$target))
    LAYER_INDEX = rep(row[1], nrow(et))
    return(cbind(et, LAYER_INDEX))
}

getLayeredEdges <- function(edgetables){
    edgetable3d = dplyr::bind_rows(edgetables)
    edgetable3d["source"] = as.character(edgetable3d$source)
    edgetable3d["target"] = as.character(edgetable3d$target)
    return(edgetable3d)
}

getLayeredNodes <- function(nodetables){
    nodetable3d = dplyr::bind_rows(nodetables)
    layeredNodes <- nodetable3d %>% rename(id = "SUID")
    layeredNodes["id"] = as.character(layeredNodes$id)
    return(layeredNodes)
}

setTransomicStyle <- function(xml, suid){
    stylename = RCy3::importVisualStyles(filename = xml)
    RCy3::setVisualStyle(stylename, network = suid)
    message(paste("Set visual style to", stylename))
}
