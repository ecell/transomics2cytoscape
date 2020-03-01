##' Import multiple KEGG pathways and integrate the pathways
##' into a cyjs file for Cy3D renderer
##'
##' @title Write cyjs file for transomics 3D visualization
##' @param pathwayID1 a KEGG pathway ID
##' @param pathwayID2 a KEGG pathway ID
##' @param pathwayID3 a KEGG pathway ID
##' @param zheight1 Z height for pathwayID1
##' @param zheight2 Z height for pathwayID2
##' @param zheight3 Z height for pathwayID3
##' @param kinase2enzyme Tsv file path for transomic interaction edges
##' @param outputcyjs Output cyjs file path
##' @author Kozo Nishida
##' @import dplyr
##' @export
##' @examples
##' library(transomics2cytoscape)
##' kinase2enzyme <- system.file("extdata", "kinase_enzyme.txt",
##'                             package = "transomics2cytoscape")
##' create3Dcyjs("rno00010", "rno00010", "rno04910", 1, 200, 400,
##'              kinase2enzyme, "transomics3D.cyjs")

create3Dcyjs <- function(pathwayID1, pathwayID2, pathwayID3, zheight1, zheight2, zheight3,
                         kinase2enzyme, outputcyjs) {
  apps=RCy3::getInstalledApps()
  #checking Apps
  if (length(grep('Cy3D,', apps)) == 0) {
    print("Cy3D is not installed yet. transomics2cytoscape installs Cy3D.")
    RCy3::installApp('Cy3D')
  }
  if (length(grep('KEGGScape,', apps)) == 0) {
    print("KEGGScape is not installed yet. transomics2cytoscape installs KEGGScape.")
    RCy3::installApp('KEGGscape')
  }
  
  writeLines(KEGGREST::keggGet(pathwayID1, option = 'kgml'), paste(pathwayID1, '.xml', sep=''))
  writeLines(KEGGREST::keggGet(pathwayID2, option = 'kgml'), paste(pathwayID2, '.xml', sep=''))
  writeLines(KEGGREST::keggGet(pathwayID3, option = 'kgml'), paste(pathwayID3, '.xml', sep=''))
  
  # layer1
  print(paste("Importing", pathwayID1))
  RCy3::importNetworkFromFile(file = paste(pathwayID1, '.xml', sep=''))
  Sys.sleep(3)
  print(paste("Reading node table", pathwayID1))
  nodetable1 = RCy3::getTableColumns()
  
  et1 = RCy3::getTableColumns(table = "edge")
  ei1 = RCy3::getEdgeInfo(et1$SUID)
  et1['source'] = unlist(lapply(ei1, function(x) x$source))
  et1['target'] = unlist(lapply(ei1, function(x) x$target))
  
  KEGG_NODE_Z = rep(zheight1, dim(nodetable1)[1])
  nodetable1 = cbind(nodetable1, KEGG_NODE_Z)
  
  # layer2
  print(paste("Importing", pathwayID2))
  RCy3::importNetworkFromFile(file = paste(pathwayID2, '.xml', sep=''))
  Sys.sleep(3)
  print(paste("Reading node table", pathwayID2))
  nodetable2 = RCy3::getTableColumns()
  
  et2 = RCy3::getTableColumns(table = "edge")
  ei2 = RCy3::getEdgeInfo(et2$SUID)
  et2['source'] = unlist(lapply(ei2, function(x) x$source))
  et2['target'] = unlist(lapply(ei2, function(x) x$target))
  
  KEGG_NODE_Z = rep(zheight2, dim(nodetable2)[1])
  nodetable2 = cbind(nodetable2, KEGG_NODE_Z)
  
  #layer3
  print(paste("Importing", pathwayID3))
  RCy3::importNetworkFromFile(file = paste(pathwayID3, '.xml', sep=''))
  Sys.sleep(3)
  print(paste("Reading node table", pathwayID3))
  nodetable3 = RCy3::getTableColumns()
  
  et3 = RCy3::getTableColumns(table = "edge")
  ei3 = RCy3::getEdgeInfo(et3$SUID)
  et3['source'] = unlist(lapply(ei3, function(x) x$source))
  et3['target'] = unlist(lapply(ei3, function(x) x$target))
  
  KEGG_NODE_Z = rep(zheight3, dim(nodetable3)[1])
  nodetable3 = cbind(nodetable3, KEGG_NODE_Z)
  
  # Combining
  print("Combining networks...")
  nodetable3d = bind_rows(nodetable1, nodetable2, nodetable3)
  nodetable3d %>% rename(id=SUID) -> nodes
  nodes['id'] = as.character(nodes$id)
  
  edges = bind_rows(et1, et2, et3)
  edges['source'] = as.character(edges$source)
  edges['target'] = as.character(edges$target)
  
  #print(kinase2enzyme)
  k2e = utils::read.table(kinase2enzyme)
  
  for(i in 1:nrow(k2e)) {
    kerow <- k2e[i,]
    source = as.character(kerow$V1)
    target = as.character(kerow$V2)
    for(j in 1:nrow(nodes)) {
      row4source = nodes[j,]
      if (source %in% row4source$KEGG_ID) {
        #print("source")
        #print(noderow$id)
        s = row4source$id
        for (k in 1:nrow(nodes)) {
          row4target = nodes[k,]
          if (target %in% row4target$KEGG_ID){
            t = row4target$id
            #print("target")
            edges %>% add_row(source = s, target = t) -> edges
            #print(c(s,t))
          }
        }
      }
    }
  }
  
  RCy3::createNetworkFromDataFrames(nodes, edges)
  
  # download.file("https://raw.githubusercontent.com/ecell/transomics2cytoscape/master/data/transomics.xml",
  #               "transomics.xml")
  stylexml <- system.file("extdata", "transomics.xml",
                               package = "transomics2cytoscape")
  #importVisualStyles(filename = "transomics.xml")
  #print(stylexml)
  RCy3::importVisualStyles(filename = stylexml)
  RCy3::setVisualStyle("transomics")
  print("Set visual style to 'transomics'")
  
  RCy3::exportNetwork(outputcyjs,'cyjs')
  print(paste('Wrote ', outputcyjs, '.cyjs in the current working directory.', sep=''))
  print("Please import it from Cytoscape Desktop and select 'Cy3D' from 'Network View Renderer' dropdown list.")
  print("Then you should see the multi layered 3D network space.")
  
  #print(nodes$KEGG_ID[[263]][3])
  #return(c(nodes, edges))
}