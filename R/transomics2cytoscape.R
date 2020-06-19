##' Import multiple KEGG pathways and integrate the pathways
##' into a cyjs file for Cy3D renderer
##'
##' @title Write cyjs file for transomics 3D visualization.
##' @param pathwayZheights A list with KEGG pathway ID and Z height in 3D
##' space.
##' @param kinase2enzyme A TSV file path between pathway network layers.
##' @param outputcyjs An output path of 3D network cyjs file.
##' @param stylexml A path of Cytoscape style file to be applied to 3D network.
##' @return This function writes a cyjs file and returns no value. 
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples
##' kinase2enzyme <- system.file("extdata", "kinase_enzyme.txt",
##'     package = "transomics2cytoscape")
##' stylexml <- system.file("extdata", "transomics.xml",
##'     package = "transomics2cytoscape")
##' create3Dcyjs(c(rno00010=1, rno00010=200, rno04910=400, rno04910=600),
##'     kinase2enzyme, stylexml, "transomics3D")

create3Dcyjs <- function(pathwayZheights, kinase2enzyme, stylexml, outputcyjs) {
    tryCatch({
        RCy3::cytoscapePing()
    }, error = function(e) {
        stop("can't connect to Cytoscape. \n
            Please check that Cytoscape is up and running.")
    })
    checkCyApps()
    pathwayIDs <- names(pathwayZheights)
    lapply(pathwayIDs, getKgml)
    suIDs <- lapply(pathwayIDs, importKgml)
    nodetables <- mapply(getNodeTableWithZheight, suIDs, pathwayZheights)
    edgetables <- lapply(suIDs, getEdgeTable)
    message("Layering the networks and adding transomics edges in 3D space...")
    layeredNodes <- getLayeredNodes(nodetables)
    transomicsEdges <- addTransomicsEdges(edgetables, kinase2enzyme,
                                            layeredNodes)
    suID <- RCy3::createNetworkFromDataFrames(layeredNodes, transomicsEdges)
    setTransomicStyle(stylexml, suID)
    exportTransomicCyjs(outputcyjs, suID)
}

checkCyApps <- function(){
    apps = RCy3::getInstalledApps()
    # checking Apps
    if (length(grep("Cy3D,", apps)) == 0) {
        message("Cy3D is not installed. transomics2cytoscape installs Cy3D.")
        RCy3::installApp("Cy3D")
    }
    if (length(grep("KEGGScape,", apps)) == 0) {
        message("KEGGScape is not installed yet.\n
        transomics2cytoscape installs KEGGScape.")
        RCy3::installApp("KEGGscape")
    }
}

getKgml <- function(pathwayID){
    writeLines(KEGGREST::keggGet(pathwayID, option = "kgml"),
                paste(pathwayID, ".xml", sep = ""))
}

importKgml <- function(pathwayID){
    message(paste("Importing", pathwayID))
    suID <- RCy3::importNetworkFromFile(file = paste(pathwayID, ".xml",
                                                        sep = ""))
    Sys.sleep(5)
    RCy3::setVisualStyle("KEGG Style")
    RCy3::fitContent()
    return(suID)
}

getNodeTableWithZheight <- function(suid, zheight){
    nodetable = RCy3::getTableColumns(table = "node", network = suid$networks)
    KEGG_NODE_Z = rep(zheight, dim(nodetable)[1])
    return(cbind(nodetable, KEGG_NODE_Z))
}

getEdgeTable <- function(suid){
    et = RCy3::getTableColumns(table = "edge", network = suid$networks)
    ei = RCy3::getEdgeInfo(et$SUID, network = suid$networks)
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

addTransomicsEdges <- function(edgetables, kinase2enzyme, layeredNodes){
    edges = bind_rows(edgetables)
    edges["source"] = as.character(edges$source)
    edges["target"] = as.character(edges$target)
    k2e = utils::read.table(kinase2enzyme)
    for (i in seq_len(nrow(k2e))) {
        kerow <- k2e[i, ]
        source = as.character(kerow$V1)
        target = as.character(kerow$V2)
        for (j in seq_len(nrow(layeredNodes))) {
            row4source = layeredNodes[j, ]
            if (source %in% unlist(row4source$KEGG_ID)) {
                s = row4source$id
                for (k in seq_len(nrow(layeredNodes))) {
                    row4target = layeredNodes[k, ]
                    if (target %in% unlist(row4target$KEGG_ID)) {
                        t = row4target$id
                        edges <- edges %>% add_row(source = s, target = t)
                    }
                }
            }
        }
    }
    return(unique(edges))
}

setTransomicStyle <- function(xml, suid){
    RCy3::importVisualStyles(filename = xml)
    RCy3::setVisualStyle("transomics", network = suid)
    message("Set visual style to 'transomics'")
}

exportTransomicCyjs <- function(outputcyjs, suid){
    RCy3::exportNetwork(outputcyjs, "cyjs", network = suid)
    cyjspath = normalizePath(paste(outputcyjs, ".cyjs", sep = ""))
    message(paste("Wrote cyjs file in ", cyjspath, sep = ""))
    message(paste("Please import", cyjspath,
    "from Cytoscape Desktop and select 'Cy3D' from 'Network View Renderer'",
    "dropdown list. Then you should see the multi layered 3D network space."))
}