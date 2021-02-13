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
    layerTable <- utils::read.table(networkLayers)
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
    # suID <- RCy3::createNetworkFromDataFrames(layeredNodes, layeredEdges)
    
    # 
    #message("Creating transomic edges to the 3D network...")
    #createTransomicsEdges(suID, transomicEdges)
    #message("Finished create3Dnetwork function!")
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
    transomicTable <- utils::read.table(transomicEdges)
    nt = RCy3::getTableColumns(table = "node", network = suid)
    et = RCy3::getTableColumns(table = "edge", network = suid)
    apply(transomicTable, 1, createTransomicEdge, nt, et)
    return(suid)
}

createTransomicEdge <- function(row, nt, et) {
    sourceLayerIndex = row[1]
    sourceTableType = row[2]
    sourceTableColumnName = row[3]
    sourceTableValue = row[4]
    targetLayerIndex = row[5]
    targetTableType = row[6]
    targetTableColumnName = row[7]
    targetTableValue = row[8]
    transomicEdgeType= row[9]
    if (sourceTableType == "node" && targetTableType == "edge"){
        createNode2Edge(nt, sourceLayerIndex, sourceTableValue,
                        sourceTableColumnName, et, targetLayerIndex,
                        targetTableValue, targetTableColumnName,
                        transomicEdgeType)
    }
    
    # if (sourceTableType == "node") {
    #     layerNt = dplyr::filter(nt, LAYER_INDEX == sourceLayerIndex)
    #     sourceNodeRows = dplyr::filter(layerNt, grepl(sourceTableValue,
    #                                     !!as.name(sourceTableColumnName)))
    #     if (targetTableType == "edge") {
    #         layerEt = dplyr::filter(et, LAYER_INDEX == targetLayerIndex)
    #         targetEdgeRows = dplyr::filter(layerEt, grepl(targetTableValue,
    #                                     !!as.name(targetTableColumnName)))
    #         if (nrow(targetEdgeRows) > 0) {
    #             ei = RCy3::getEdgeInfo(targetEdgeRows["SUID"])
    #             #centerNodes = lapply(ei, createNodeForEdge)
    #             reactionSourceNodes = lapply(ei, getEdgeSourceSUID)
    #             
    #             for (i in seq_len(nrow(sourceNodeRows))){
    #                 sourceSUID = sourceNodeRows[i, 1]
    #                 for (j in seq_len(length(reactionSourceNodes))){
    #                     targetSUID = reactionSourceNodes[[j]]
    #                     RCy3::addCyEdges(c(sourceSUID, targetSUID),
    #                         edgeType=as.character(transomicEdgeType))
    #                 }
    #             }
    #         }
    #     }
    # }
}

createNode2Edge <- function(nt, sourceLayerIndex, sourceTableValue,
                            sourceTableColumnName, et, targetLayerIndex,
                            targetTableValue, targetTableColumnName,
                            transomicEdgeType){
    layerNt = dplyr::filter(nt, LAYER_INDEX == sourceLayerIndex)
    sourceNodeRows = dplyr::filter(layerNt, grepl(sourceTableValue,
                                    !!as.name(sourceTableColumnName)))
    layerEt = dplyr::filter(et, LAYER_INDEX == targetLayerIndex)
    targetEdgeRows = dplyr::filter(layerEt, grepl(targetTableValue,
                                    !!as.name(targetTableColumnName)))
    if (nrow(targetEdgeRows) > 0) {
        ei = RCy3::getEdgeInfo(targetEdgeRows["SUID"])
        #centerNodes = lapply(ei, createNodeForEdge)
        reactionSourceNodes = lapply(ei, getEdgeSourceSUID)
        for (i in seq_len(nrow(sourceNodeRows))){
            sourceSUID = sourceNodeRows[i, 1]
            for (j in seq_len(length(reactionSourceNodes))){
                targetSUID = reactionSourceNodes[[j]]
                RCy3::addCyEdges(c(sourceSUID, targetSUID),
                                 edgeType=as.character(transomicEdgeType))
            }
        }
    }
}

getEdgeSourceSUID <- function(edgeInfo){
    sourceNodeSUID = edgeInfo$source
    return(sourceNodeSUID)
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
    transomicTable = utils::read.table(tsvFilePath)
    ecVec = transomicTable[ , columnIndex]
    ecVec = unique(ecVec)
    ec2rea = KEGGREST::keggLink("reaction", ecVec)
    dflist = apply(transomicTable, 1, ecRow2reaRows, columnIndex, ec2rea)
    df = dplyr::bind_rows(dflist)
    df = dplyr::select(df, -one_of("rows"))
    lastColumn = df[ , ncol(df)]
    originalColumn = df[ , columnIndex]
    df[ , columnIndex] = lastColumn
    df[ , ncol(df)] = originalColumn
    utils::write.table(df, file=outputFilename, quote=FALSE, sep='\t',
                        col.names = F, row.names = F)
}

ecRow2reaRows <- function(row, columnIndex, ec2rea) {
    ec = row[columnIndex]
    rea = ec2rea[names(ec2rea) == ec]
    rows = do.call("rbind", replicate(length(rea), row, simplify = FALSE))
    return(as.data.frame(cbind(rows, as.vector(rea))))
}

# createNodeForEdge <- function(edgeInfo){
#     sourceNodeSUID = edgeInfo$source
#     targetNodeSUID = edgeInfo$target
#     theName = paste(as.character(sourceNodeSUID), as.character(targetNodeSUID))
#     newNodeInfo = RCy3::addCyNodes(theName, skip.duplicate.names = TRUE)
#     if (length(newNodeInfo) > 0) {
#         newNodeSUID = newNodeInfo[[1]]$SUID
#         # RCy3::addCyEdges(list(c(sourceNodeSUID, newNodeSUID),
#         #                       c(newNodeSUID, targetNodeSUID)))
#         return(c(newNodeSUID, sourceNodeSUID, targetNodeSUID))    
#     }
# }

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
        networkSUID = res$networks
        return(networkSUID)
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
    RCy3::importVisualStyles(filename = xml)
    RCy3::setVisualStyle("transomics", network = suid)
    message("Set visual style to 'transomics'")
}
