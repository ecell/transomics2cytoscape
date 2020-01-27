library(RCy3)
library(dplyr)

installApp('KEGGscape')

create3Dnetfile <- function(kgml1, kgml2, kgml3, zheight1, zheight2, zheight3) {
  importNetworkFromFile(file = kgml1) 
  importNetworkFromFile(file = kgml2)
  importNetworkFromFile(file = kgml3)
  
  nodetable1 = getTableColumns('node', network = getNetworkList()[1])
  KEGG_NODE_Z = rep(zheight1, dim(nodetable1)[1])
  nodetable1 = cbind(nodetable1, KEGG_NODE_Z)
  
  nodetable2 = getTableColumns('node', network = getNetworkList()[2])
  KEGG_NODE_Z = rep(zheight2, dim(nodetable2)[1])
  nodetable2 = cbind(nodetable2, KEGG_NODE_Z)  
  
  nodetable3 = getTableColumns('node', network = getNetworkList()[3])
  KEGG_NODE_Z = rep(zheight3, dim(nodetable3)[1])
  nodetable3 = cbind(nodetable3, KEGG_NODE_Z)
  
  nodetable3d = bind_rows(nodetable1, nodetable2, nodetable3)
  nodetable3d %>% rename(id=SUID) -> nodes
  nodes['id'] = as.character(nodes$id)
  createNetworkFromDataFrames(nodes)
  exportNetwork('./transomics','cyjs')
}

## example
#create3Dnetfile('rno00010.xml', 'rno00010.xml', 'rno04910.xml', 1, 500, 1000)
