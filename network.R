# inja man mikham ye network bekesham ke vertex ha bimari ha va mir ha bashan
# edge bine bimari X va mir Y hast agar va tanha agar , mir Y  dif.ex cancer X va tamami normal haye dge manfi bashe 
# va significant ham bashe 
# dar vaghe edge ha az anwer tabe topMirs_AllNomrals miad ke dakhel findMir.R hast.


tis_names = c("BLCA" , "BRCA" , "ESCA" , "HNSC" , "KICH" , "KIRC" , "KIRP" , "LIHC" , "LUAD" , "LUSC" , "PRAD" , "STAD" , "THCA" , "UCEC" )

#library(RCytoscape)
g <- new ("graphNEL", edgemode="directed")
cw <- new.CytoscapeWindow ("my_net", graph=g)



# preprocess 
setOfMirs = c()
for ( i in 1:14 ) {
  tmpName = tis_names[i];
  #list = topMirs_AllNomrals(name=tmpName,num=10)
  list = topMirs_AllNomrals(name=tmpName)
  setOfMirs = unique(c(setOfMirs,list))
}


#------------------ node haye graph --------------------------

# dasteE az nod ha bimari ha hastan

for ( i in 1:14 ) {
  name_ = tis_names[i]
  g <- graph::addNode (as.character(name_), g)
  setNodeShapeDirect(cw,name_,new.shapes="Diamond")
  setNodeSizeDirect(cw,name_,80)
  setNodeLabelColorDirect(cw,name_,"#000000")
  setNodeColorDirect(cw,name_,"#06B2A5")
}

# daste digar mir ha 

for ( i in 1:length(setOfMirs) ) {
  name_ = as.character(setOfMirs[i])
  g <- graph::addNode (as.character(name_), g)
  setNodeSizeDirect(cw,name_,25)
  setNodeLabelColorDirect(cw,name_,"#000000")
  setNodeColorDirect(cw,name_,"#67204A")
}
cw = setGraph (cw, g)
displayGraph (cw)

layoutNetwork (cw, layout.name="grid")
redraw (cw)





#-------------------------- add edge --------------------------

for ( i in 1:14 ) {
  tmpName =  tis_names[i];
  list = topMirs_AllNomrals(name=tmpName)
  for ( j in 1:length(list) ) {
    g <- graph::addEdge (as.character(tis_names[i]), as.character(list[j]), g)
  }
  
}




#-------------------------- size node ha ---------------------- 
for ( i in 1:length(setOfMirs) ) {
  name_ = as.character(setOfMirs[i])
  
  if (as.numeric(degree(g,name_)[1]) == 1) {
    setNodeSizeDirect(cw,name_,60)
  } 
}

cw = setGraph (cw, g)
displayGraph (cw)    # cw's graph is sent to Cytoscape
redraw (cw)






