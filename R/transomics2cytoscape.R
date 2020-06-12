##' Import multiple KEGG pathways and integrate the pathways
##' into a cyjs file for Cy3D renderer
##'
##' @title Write cyjs file for transomics 3D visualization
##' @param pathwayList a list with KEGG pathway ID and the Z height in 3D space
##' @param kinase2enzyme Tsv file path for transomic interaction edges
##' @param outputcyjs Output cyjs file path
##' @return This function writes a cyjs file and returns no value. 
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples
##' kinase2enzyme <- system.file('extdata', 'kinase_enzyme.txt',
##'                             package = 'transomics2cytoscape')
##' create3Dcyjs(c(rno00010=1, rno00010=200, rno04910=400),
##'              kinase2enzyme, 'transomics3D')

create3Dcyjs <- function(pathwayList, kinase2enzyme, outputcyjs) {
  checkCyApps()
  
  pathwayIDs <- names(pathwayList)
  lapply(pathwayIDs, getKgml)
  suIDs <- lapply(pathwayIDs, importKgml)
  nodetables <- mapply(getNodeTableWithZheight, suIDs, pathwayList)
  edgetables <- lapply(suIDs, getEdgeTable)
  
  # Combining
  print("Combining networks...")
  print(length(nodetables))
  nodetable3d = bind_rows(nodetables)
  nodes <- nodetable3d %>% rename(id = "SUID")
  nodes["id"] = as.character(nodes$id)
  
  edges = bind_rows(edgetables)
  edges["source"] = as.character(edges$source)
  edges["target"] = as.character(edges$target)
  
  # print(kinase2enzyme)
  k2e = utils::read.table(kinase2enzyme)
  
  # return(list(k2e, nodes, edges))
  
  for (i in 1:nrow(k2e)) {
    kerow <- k2e[i, ]
    source = as.character(kerow$V1)
    target = as.character(kerow$V2)
    for (j in 1:nrow(nodes)) {
      row4source = nodes[j, ]
      if (source %in% unlist(row4source$KEGG_ID)) {
        # print('source') print(noderow$id)
        s = row4source$id
        for (k in 1:nrow(nodes)) {
          row4target = nodes[k, ]
          if (target %in% unlist(row4target$KEGG_ID)) {
            t = row4target$id
            edges <- edges %>% add_row(source = s, target = t)
          }
        }
      }
    }
  }
  
  RCy3::createNetworkFromDataFrames(nodes, unique(edges))
  
  stylexml <- system.file("extdata", "transomics.xml",
                          package = "transomics2cytoscape")
  RCy3::importVisualStyles(filename = stylexml)
  RCy3::setVisualStyle("transomics")
  print("Set visual style to 'transomics'")
  
  RCy3::exportNetwork(outputcyjs, "cyjs")
  print(paste("Wrote ", outputcyjs,
              ".cyjs in the current working directory.", sep = ""))
  print("Please import it from Cytoscape Desktop and select 'Cy3D' from\n
        'Network View Renderer' dropdown list.")
  print("Then you should see the multi layered 3D network space.")
  
  # print(nodes$KEGG_ID[[263]][3]) return(c(nodes, edges))
}

checkCyApps <- function(){
  apps = RCy3::getInstalledApps()
  # checking Apps
  if (length(grep("Cy3D,", apps)) == 0) {
    print("Cy3D is not installed yet. transomics2cytoscape installs Cy3D.")
    RCy3::installApp("Cy3D")
  }
  if (length(grep("KEGGScape,", apps)) == 0) {
    print("KEGGScape is not installed yet.\n
    transomics2cytoscape installs KEGGScape.")
    RCy3::installApp("KEGGscape")
  }
}

getKgml <- function(pathwayID){
  writeLines(KEGGREST::keggGet(pathwayID, option = "kgml"),
             paste(pathwayID, ".xml", sep = ""))
}

importKgml <- function(pathwayID){
  print(paste("Importing", pathwayID))
  suID <- RCy3::importNetworkFromFile(file = paste(pathwayID, ".xml", sep = ""))
  Sys.sleep(3)
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