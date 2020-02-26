installApp('KEGGscape')
installApp('Cy3D')

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
##' @param stylexml Cytoscape style XML file path for 3D visualization
##' @param kinase2enzyme Tsv file path for transomic interaction edges
##' @param outputcyjs Output cyjs file path
##' @author Kozo Nishida
##' @export
##' @examples
##' create3Dcyjs()

create3Dcyjs <- function(pathwayID1, pathwayID2, pathwayID3, zheight1, zheight2, zheight3,
                         stylexml, kinase2enzyme, outputcyjs) {
  
  writeLines(keggGet(pathwayID1, option = 'kgml'), paste(pathwayID1, '.xml', sep=''))
  writeLines(keggGet(pathwayID2, option = 'kgml'), paste(pathwayID2, '.xml', sep=''))
  writeLines(keggGet(pathwayID3, option = 'kgml'), paste(pathwayID3, '.xml', sep=''))
  
  # layer1
  print(paste("Importing", pathwayID1))
  importNetworkFromFile(file = paste(pathwayID1, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID1))
  nodetable1 = getTableColumns()
  
  et1 = getTableColumns(table = "edge")
  ei1 = getEdgeInfo(et1$SUID)
  et1['source'] = unlist(lapply(ei1, function(x) x$source))
  et1['target'] = unlist(lapply(ei1, function(x) x$target))
  
  KEGG_NODE_Z = rep(zheight1, dim(nodetable1)[1])
  nodetable1 = cbind(nodetable1, KEGG_NODE_Z)
  
  # layer2
  print(paste("Importing", pathwayID2))
  importNetworkFromFile(file = paste(pathwayID2, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID2))
  nodetable2 = getTableColumns()
  
  et2 = getTableColumns(table = "edge")
  ei2 = getEdgeInfo(et2$SUID)
  et2['source'] = unlist(lapply(ei2, function(x) x$source))
  et2['target'] = unlist(lapply(ei2, function(x) x$target))
  
  KEGG_NODE_Z = rep(zheight2, dim(nodetable2)[1])
  nodetable2 = cbind(nodetable2, KEGG_NODE_Z)
  
  #layer3
  print(paste("Importing", pathwayID3))
  importNetworkFromFile(file = paste(pathwayID3, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID3))
  nodetable3 = getTableColumns()
  
  et3 = getTableColumns(table = "edge")
  ei3 = getEdgeInfo(et3$SUID)
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
  
  k2e = read.table(kinase_enzyme)
  
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
  
  createNetworkFromDataFrames(nodes, edges)
  
  download.file("https://raw.githubusercontent.com/ecell/transomics2cytoscape/master/transomics.xml", "transomics.xml")
  importVisualStyles(filename = "transomics.xml")
  setVisualStyle("transomics")
  print("Set visual style to 'transomics'")
  
  exportNetwork(output,'cyjs')
  print(paste('Wrote ', output, '.cyjs in the current working directory.', sep=''))
  print("Please import it from Cytoscape Desktop and select 'Cy3D' from 'Network View Renderer' dropdown list.")
  print("Then you should see the multi layered 3D network space.")
  
  #print(nodes$KEGG_ID[[263]][3])
  #return(c(nodes, edges))
}