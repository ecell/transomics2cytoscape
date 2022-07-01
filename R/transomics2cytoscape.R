##' Import multiple KEGG pathways and integrate the pathways
##' into Cy3D renderer
##'
##' @title Create 3D network view for transomics visualization.
##' @param networkDataDir Path of a directory to put the network files
##' of the second column of networkLayers TSV.
##' @param networkLayers Path of a TSV file with the 4 columns (layer index,
##' the network file name in networkDataDir, Z-height of the network,
##' whether to interact not only with the nodes of each network layer but also
##' with the edges).
##' @param stylexml Path of a XML file for Cytoscape style
##' @return A SUID of the 3D network.
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' networkDataDir <- tempfile(); dir.create(networkDataDir)
##' networkLayers <- system.file("extdata/usecase1", "yugi2014.tsv",
##'     package = "transomics2cytoscape")
##' stylexml <- system.file("extdata/usecase1", "yugi2014.xml",
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
    layerTable = utils::read.table(networkLayers, sep="\t")
    networkSUID = apply(layerTable, 1, importLayer)
    layerTable = cbind(layerTable, networkSUID)
    nodetables <- apply(layerTable, 1, getNodeTableWithLayerinfo)
    layeredNodes <- getLayeredNodes(nodetables)
    edgetables <- apply(layerTable, 1, getEdgeTableWithLayerinfo)
    message("I'm in checkpoint 1 ...")
    layeredEdges <- getLayeredEdges(edgetables)
    message("I'm in checkpoint 2 ...")
    return(c(layeredNodes, layeredEdges))
    #RCy3::commandsPOST('cy3d set renderer')
    #suID <- RCy3::createNetworkFromDataFrames(layeredNodes, layeredEdges)
    #setTransomicStyle(stylexml, suID)
    #return(suID)
}

##' Create Trans-Omic edges between layers of the network
##'
##' @title Create Trans-Omic edges between layers of the network.
##' @param suid A SUID of Cytoscape network
##' @param transomicEdges Path of a TSV file with the 7 columns
##' (layer index of the source node,
##' the column name for which you want to find the attribute value of the
##' source node,
##' the attribute value of the source node should have,
##' layer index of a target node
##' name or KEGG object ID that the source node should have,
##' layer index of the target node,
##' the column name for which you want to find the attribute value of the
##' target node,
##' the attribute value of the target node should have,
##' name or KEGG object ID that the target node should have,
##' interaction type).
##' @return A SUID of the 3D network.
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' layer1to2 <- system.file("extdata/usecase1", "k2e.tsv",
##'     package = "transomics2cytoscape")
##' suid <- createTransomicEdges(suid, layer1to2)
##' }

createTransomicEdges <- function(suid, transomicEdges) {
    transomicTable <- utils::read.table(transomicEdges, sep="\t")
    nt = RCy3::getTableColumns(table = "node", network = suid)
    #et = RCy3::getTableColumns(table = "edge", network = suid)
    #addedEdges = list()
    edgetbl = tibble(source=numeric(), target=numeric(),
                    interaction=character())
    for (i in seq_len(nrow(transomicTable))) {
        row = transomicTable[i,]
        #edgetbl = createTransomicEdge(row, nt, et, edgetbl, suid)
        edgetbl = createTransomicEdge(row, nt, edgetbl, suid)
    }
    #return(edgetbl)
    body = purrr::transpose(dplyr::distinct(edgetbl))
    RCy3::cyrestPOST(
        paste("networks", suid, "edges", sep = "/"),
        body = body
    )
    return(suid)
}

createTransomicEdge <- function(row, nt, edgetbl, suid) {
    sourceLayerIndex = as.character(row[1])
    #sourceTableType = as.character(row[2])
    sourceTableColumnName = as.character(row[2])
    sourceTableValue = as.character(row[3])
    targetLayerIndex = as.character(row[4])
    #targetTableType = as.character(row[6])
    targetTableColumnName = as.character(row[5])
    targetTableValue = as.character(row[6])
    transomicEdgeType= as.character(row[7])
    # if (sourceTableType == "node" && targetTableType == "edge"){
    #     addedEdges = createNode2Edge(nt, sourceLayerIndex, sourceTableValue,
    #                     sourceTableColumnName, et, targetLayerIndex,
    #                     targetTableValue, targetTableColumnName,
    #                     transomicEdgeType, addedEdges, suid)
    #} else if (sourceTableType == "node" && targetTableType == "node"){
    edgetbl = createNode2Node(nt, sourceLayerIndex, sourceTableValue,
                        sourceTableColumnName, targetLayerIndex,
                        targetTableValue, targetTableColumnName,
                        transomicEdgeType, edgetbl)
    # } else if (sourceTableType == "edge" && targetTableType == "edge"){
    #     addedEdges = createEdge2Edge(et, sourceLayerIndex, sourceTableValue,
    #                     sourceTableColumnName, targetLayerIndex,
    #                     targetTableValue, targetTableColumnName,
    #                     transomicEdgeType, addedEdges, suid)
    # }
    return(edgetbl)
}

createNode2Node <- function(nt, sourceLayerIndex, sourceTableValue,
                            sourceTableColumnName, targetLayerIndex,
                            targetTableValue, targetTableColumnName,
                            transomicEdgeType, edgetbl) {
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
                #if (!(list(c(sourceSUID, targetSUID)) %in% edgetbl)) {
                edgetbl = edgetbl %>% tibble::add_row(source=sourceSUID,
                            target=targetSUID,
                            interaction=as.character(transomicEdgeType))
                    # addedEdges = append(addedEdges,
                    #                     list(c(sourceSUID, targetSUID)))
                    # RCy3::addCyEdges(c(sourceSUID, targetSUID),
                    #                 edgeType=as.character(transomicEdgeType))
                #}
            }
        }
    }
    return(edgetbl)
}

##' Convert a EC number column to a KEGG reaction ID column
##'
##' @title Convert KEGG enzyme IDs to KEGG reaction IDs.
##' @param tsvFilePath Path of a TSV file with column of EC number
##' @param columnIndex Index number of the column with the EC number you want
##' to convert
##' @param outputFilename The output file name
##' @return None
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples \dontrun{
##' layer3to2 <- system.file("extdata/usecase1", "allosteric_ecnumber.tsv",
##'     package = "transomics2cytoscape")
##' ec2reaction(layer3to2, 6, "allosteric_ec2rea.tsv")
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

importLayer <- function(row){
    fileExtension <- tools::file_ext(row[2])
    if (fileExtension %in% c("sif", "gml", "xgmml", "xml")){
        res <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                        row[2], sep=""))
        Sys.sleep(3)
        nodePositions <- RCy3::getNodePosition(network = res$networks)
        names(nodePositions)[1] <- "x_location"
        names(nodePositions)[2] <- "y_location"
        RCy3::loadTableData(nodePositions, network = res$networks)
        if (row[4] == "true") {
            message("Creating a node at each midpoint of the edge...")
            createMidnodes(res$networks)
        }
        return(res$networks)
    } else {
        getKgml(row[2])
        kgml <- paste(row[2], ".xml", sep = "")
        message("Importing ", kgml)
        res <- RCy3::importNetworkFromFile(file = paste(getwd(), "/",
                                                        kgml, sep=""))
        Sys.sleep(3)
        RCy3::setVisualStyle("KEGG Style")
        xy = RCy3::getTableColumns(table = 'node',
                            columns = c('KEGG_NODE_X', 'KEGG_NODE_Y'))
        colnames(xy) = c('x_location', 'y_location')
        RCy3::loadTableData(xy, table.key.column = "SUID")
        RCy3::updateStyleMapping(RCy3::getCurrentStyle(),
                RCy3::mapVisualProperty("NODE_X_LOCATION", "x_location", "p"))
        RCy3::updateStyleMapping(RCy3::getCurrentStyle(),
                RCy3::mapVisualProperty("NODE_Y_LOCATION", "y_location", "p"))
        RCy3::fitContent()
        if (row[4] == "true") {
            message("Creating a node at each midpoint of the edge...")
            createMidnodes(res$networks)
        }
        return(res$networks)
    }
}

createMidnodes <- function(networkSUID){
    #nt = RCy3::getTableColumns(table = 'node', network = networkSUID)
    et = RCy3::getTableColumns(table = 'edge', network = networkSUID)
    message("getting node positions...")
    # nodePosition = RCy3::getNodePosition(unlist(nt[['SUID']]),
    #                                    network = networkSUID)
    xyloc = RCy3::getTableColumns(table = 'node',
        columns = c('SUID', 'x_location', 'y_location'), network = networkSUID)
    message("getting edge information...")
    einfo = RCy3::getEdgeInfo(unlist(et[['SUID']]), network = networkSUID)
    einfodf = dplyr::bind_rows(einfo)
    midxydf = dplyr::bind_rows(lapply(einfo, getMidLoc, xyloc))
    midNodeInfo = dplyr::bind_cols(einfodf, midxydf)
    midNodeInfo = dplyr::rename(midNodeInfo, orig_edge_SUID=SUID)
    message("creating midpoint nodes...")
    midNodes = RCy3::addCyNodes(as.character(et[['SUID']]))

    RCy3::loadTableData(midNodeInfo, data.key.column = "orig_edge_SUID",
                        table = "node", table.key.column = "name")
    # RCy3::loadTableData(xyloc, data.key.column = "SUID",
    #                     table = "node", table.key.column = "SUID")
    # message("updating style...")
    RCy3::updateStyleMapping(RCy3::getCurrentStyle(),
                RCy3::mapVisualProperty("NODE_X_LOCATION", "x_location", "p"))
    RCy3::updateStyleMapping(RCy3::getCurrentStyle(),
                RCy3::mapVisualProperty("NODE_Y_LOCATION", "y_location", "p"))
    message("deleting edges...")
    RCy3::invertEdgeSelection()
    RCy3::deleteSelectedEdges()
    message("start collecting node IDs to connect...")
    midNodesName2suid = dplyr::bind_rows(midNodes)
    midNodeInfo = dplyr::bind_cols(midNodeInfo, midNodesName2suid['SUID'])
    src2mid = midNodeInfo[c('source', 'SUID')]
    src2mid = dplyr::rename(src2mid, target = SUID)
    src2mid = src2mid %>% purrr::transpose()
    mid2tgt = midNodeInfo[c('SUID', 'target')]
    mid2tgt = dplyr::rename(mid2tgt, source = SUID)
    mid2tgt = mid2tgt %>% purrr::transpose()
    res = RCy3::cyrestPOST(paste("networks", networkSUID, "edges", sep = "/"),
            body=src2mid)
    #res = dplyr::bind_rows(res)
    #res = dplyr::bind_cols(res['SUID'], midNodeInfo)
    res = RCy3::cyrestPOST(paste("networks", networkSUID, "edges", sep = "/"),
            body=mid2tgt)
    #message("creating edges for midpoionts...")
    #saveRDS(edgelist, file = "midpoint_edgesrctgt.rds")
    #RCy3::addCyEdges(edgelist)
}

getMidLoc = function(elem, xyloc){
    sourcexy=xyloc[rownames(xyloc) == elem$source, ]
    targetxy=xyloc[rownames(xyloc) == elem$target, ]
    midx_location=(as.numeric(sourcexy['x_location']) +
                    as.numeric(targetxy['x_location'])) / 2
    midy_location=(as.numeric(sourcexy['y_location']) +
                    as.numeric(targetxy['y_location'])) / 2
    #midy_location=(sourcexy['y_location']+targetxy['y_location'])/2
    return(list(x_location=unlist(midx_location),
                y_location=unlist(midy_location)))
    # return(list(edge_SUID=elem$SUID, x_location=unlist(midx_location),
    #             y_location=unlist(midy_location), source=elem$source,
    #             target=elem$target))
}

connectMidnodeEdges = function(elem, midNodeInfo){
    midnodeSuid = elem$SUID
    midNodeRow = midNodeInfo[midNodeInfo['orig_edge_SUID']==elem$SUID, ]
    sourceSuid = midNodeRow['source']
    targetSuid = midNodeRow['target']
    return(c(sourceSuid, midnodeSuid))
    # return(
    #     c(setNames(c(sourceSuid, midnodeSuid), c("source", "target")),
    #       directed = FALSE, interaction = 'interacts with'))
    #return(list(c(sourceSuid,midnodeSuid),c(midnodeSuid,targetSuid)))
    #RCy3::addCyEdges(list(c(sourceSuid,midnodeSuid),c(midnodeSuid,targetSuid)))
}

getKgml <- function(pathwayID){
    writeLines(KEGGREST::keggGet(pathwayID, option = "kgml"),
                paste(pathwayID, ".xml", sep = ""))
}

getNodeTableWithLayerinfo <- function(row){
    nt = RCy3::getTableColumns(table = "node", network = as.numeric(row[5]))
    #KEGG_NODE_Z = as.character(rep(row[3], nrow(nt)))
    z_location = as.character(rep(row[3], nrow(nt)))
    LAYER_INDEX = rep(row[1], nrow(nt))
    return(cbind(nt, z_location, LAYER_INDEX))
}

getEdgeTableWithLayerinfo <- function(row){
    et = RCy3::getTableColumns(table = "edge", network = as.numeric(row[5]))
    message("Getting edge info. This function is kinda slow...")
    ei = RCy3::getEdgeInfo(et$SUID, as.numeric(row[5]))
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
