library(RCy3)
library(dplyr)
library(KEGGREST)

installApp('KEGGscape')
installApp('Cy3D')

create3Dnetfile <- function(pathwayID1, pathwayID2, pathwayID3, zheight1, zheight2, zheight3, output) {
  writeLines(keggGet(pathwayID1, option = 'kgml'), paste(pathwayID1, '.xml', sep=''))
  writeLines(keggGet(pathwayID2, option = 'kgml'), paste(pathwayID2, '.xml', sep=''))
  writeLines(keggGet(pathwayID3, option = 'kgml'), paste(pathwayID3, '.xml', sep=''))
  
  print(paste("Importing", pathwayID1))
  importNetworkFromFile(file = paste(pathwayID1, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID1))
  nodetable1 = getTableColumns()
  KEGG_NODE_Z = rep(zheight1, dim(nodetable1)[1])
  nodetable1 = cbind(nodetable1, KEGG_NODE_Z)
  
  print(paste("Importing", pathwayID2))
  importNetworkFromFile(file = paste(pathwayID2, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID2))
  nodetable2 = getTableColumns()
  KEGG_NODE_Z = rep(zheight2, dim(nodetable2)[1])
  nodetable2 = cbind(nodetable2, KEGG_NODE_Z)
  
  print(paste("Importing", pathwayID3))
  importNetworkFromFile(file = paste(pathwayID3, '.xml', sep=''))
  Sys.sleep(5)
  print(paste("Reading node table", pathwayID3))
  nodetable3 = getTableColumns()
  KEGG_NODE_Z = rep(zheight3, dim(nodetable3)[1])
  nodetable3 = cbind(nodetable3, KEGG_NODE_Z)
  
  print("Combining networks...")
  nodetable3d = bind_rows(nodetable1, nodetable2, nodetable3)
  nodetable3d %>% rename(id=SUID) -> nodes
  nodes['id'] = as.character(nodes$id)
  createNetworkFromDataFrames(nodes)
  
  download.file("https://raw.githubusercontent.com/ecell/cytoscape-styles/master/xml/transomics.xml", "transomics.xml")
  importVisualStyles(filename = "transomics.xml")
  setVisualStyle("transomics")
  print("Set visual style to 'transomics'")
  
  exportNetwork(output,'cyjs')
  print(paste('Wrote ', output, '.cyjs in the current working directory.', sep=''))
  print("Please import it from Cytoscape Desktop and select 'Cy3D' from 'Network View Renderer' dropdown list.")
  print("Then you should see the multi layered 3D network space.")
}

## example
#create3Dnetfile('rno00010', 'rno00010', 'rno04910', 1, 500, 1000, 'output')
