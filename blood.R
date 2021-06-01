library(GEOquery)




gse <- getGEO("GSE31309",GSEMatrix=FALSE) # download serie
GSMList(gse)[[1]] -> gsm1 # seprate samples
head(Table(gsm1)) # table for samples but also we have meta data for every of them.


makeTableFromBloodSeris <- function( GSE , nids , nide, tids , tide ) {
  gse = getGEO(GSE , GSEMatrix=FALSE)
  platform = GPLList(gse)
  annt = Table(platform[[1]])
  samples = GSMList(gse)
  samplest = lapply(samples , Table)
  table = do.call(cbind,samplest)
  table$ID = table[,1]
  table = table[,-grep("REF" , colnames(table))]
  
  fixAnnotation( table, annt ) -> table
  
  coln = sapply( as.numeric(nids):as.numeric(nide) , function(x) { t = colnames(table)[x]; t = strsplit2(t,split="\\.")[1] ; paste(t,"Normal" , sep=".") } )
  colt = sapply( as.numeric(tids):as.numeric(tide) , function(x) { t = colnames(table)[x]; t = strsplit2(t,split="\\.")[1] ; paste(t,"Tumor" , sep=".") } )
  colnames(table) <- c(coln,colt)
  table
}
fixAnnotation <- function( table , annt ) {
  merge(table , annt , by.x = "ID" , by.y = "ID") -> table
  table$ID = NULL
  table$ID_description = NULL
  table$SEQUENCE = NULL
  row.names(table) <- table$miRNA_ID
  table$miRNA_ID = NULL
  table
}




########## analys data 




PCA_oneType <-function( myData , name , index  ) {
  pc = prcomp(t(myData))
  pcx = as.data.frame(pc$x)
  Sample = substr(row.names(pcx) , index, index)
  
  pdf(paste( "~/Desktop/Saarland/results/blood/pca" , paste (name , ".pdf", sep="") , sep= "_"), width=12, height=8)
  print(ggplot(data = pcx, aes(x = PC1, y = PC2, color=Sample)) +
          geom_point()+ 
          theme(panel.background = element_rect(fill = 'white', colour = 'red'),
                panel.border = element_rect(fill = NA, colour = "black", size=2),
                plot.background = element_rect(fill = '#FFE4E1'),
                strip.text.x = element_text(size = 25, colour = "black"),
                strip.background = element_rect(fill = "#FFE4E1"),
                axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
                axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
                legend.background = element_rect(fill = '#FFE4E1',colour=NA), 
                legend.key =        element_rect(fill =NA, colour = "black", size = 0.25),
                legend.key.size =   unit(1.5, "lines"),
                legend.text =       element_text(size = 20 * 1),
                legend.title =      element_text(size = 20 * 1),
                legend.position =   "right"
                
          )
        ) 
  dev.off()

}


############################# calculate DE for bloods
tt = makeTableFromBloodSeris("GSE31309",1,57,58,105)
tt$ID = rownames(tt)
s2 = read.csv2("~/Desktop/Saarland/conv/converted_mirnas.csv", sep = "\t" , header = T)
tt = merge( tt , s2 , by.x = "ID" , by.y = "v15")
tt$ID = NULL
tt = tt[!duplicated(tt$v21),]
row.names(tt) <- tt$v21
tt$v21 = NULL

data = tt

ex = data
gr = substr(colnames(ex),1,1)
con = c("N-T")
answer = dif.ex(ex,gr,con)
ans = answer[order(answer$logFC,decreasing= F) , ]
ans = ans[ which ( ans$logFC < -1.5) , ]
ans = ans[ which ( ans$adj.P.Val < 0.01) , ]
ans = ans[order(ans$logFC,decreasing=F),]


##################################################### PCA for 54 mirs DE from BRCA TCGA
tt = makeTableFromBloodSeris("GSE31309",1,57,58,105)
tt = normalizeQuantiles(tt)
tt$ID = rownames(tt)
s2 = read.csv2("~/Desktop/Saarland/conv/converted_mirnas.csv", sep = "\t" , header = T)
tt = merge( tt , s2 , by.x = "ID" , by.y = "v15")
tt$ID = NULL
tt = tt[!duplicated(tt$v21),]
row.names(tt) <- tt$v21
tt$v21 = NULL
#rr is rownames DE for BRCA 
tt = tt[rr,]
tt = tt[-which( is.na(tt[,1]) == TRUE) , ]
myData = tt


########################################### DE baraye jofteshun hamzaman,
tt = makeTableFromBloodSeris("GSE31309",1,57,58,105)
tt$ID = rownames(tt)
s2 = read.csv2("~/Desktop/Saarland/conv/converted_mirnas.csv", sep = "\t" , header = T)
tt = merge( tt , s2 , by.x = "ID" , by.y = "v15")
tt$ID = NULL
tt = tt[!duplicated(tt$v21),]
row.names(tt) <- tt$v21
tt$v21 = NULL
bdata = tt
bdata$ID = row.names(bdata)

tdata = read.delim("~/Desktop/Saarland/data/BRCA.txt")
tdata = normalizeQuantiles(tdata)
tdata$ID = row.names(tdata)
mdata = merge (tdata , bdata , by.x = "ID" , by.y = "ID", all = FALSE )
rownames(mdata) <- mdata$ID
mdata$ID = NULL


ex = mdata
gr = c(substr(colnames(ex[1:865]),1,1) , paste0("B",substr(colnames(ex[866:970]),11,11)) )
con = c("N-T" , "BN-BT")
answer = dif.ex(ex,gr,con)
ans = answer
#ans = ans[ which ( ans$logFC < 0) , ]
ans = ans[ which ( ans$adj.P.Val < 0.01) , ]
ans = ans[ which ( abs( ans$BN.BT) > 0.4) , ]
ans = ans[ which ( abs( ans$N.T)  > 0.6) , ]
ans = ans[ which ( ans$N.T * ans$BN.BT > 0) ,]
ans = ans[order(ans$N.T,decreasing=F),]

### Averages : 
row.names(ans) -> rn

bdata2 = bdata[rn,]
bdata2$ID = NULL
aveBloodNoraml = as.data.frame(t(t( apply(bdata2[,1:57],1,mean))))
aveBloodTumor = as.data.frame(t(t( apply(bdata2[,58:105],1,mean) ) ) ) 

tdata2 = tdata[rn,]
tdata2$ID = NULL
aveTissueNormal = as.data.frame(t(t( apply(tdata2[,779:865],1,mean) )))
aveTissueTumor = as.data.frame(t(t( apply(tdata2[,1:778],1,mean) )))

myMean = as.data.frame(t(t(myMean)))
myMean = cbind(myMean,name=row.names(myMean))

###########
ans = answer
ansUU = ans[which ( ans$BN.BT < -0.2) , ]
ansUU = ansUU[ which( ansUU$N.T < -0.2 ) , ]

ansUD = ans[ which( ans$BN.BT < -0.2) ,]
ansUD = ansUD[ which( ansUD$N.T > 0.2) ,]

ansDU = ans[ which( ans$BN.BT > 0.2) ,]
ansDU = ansDU[ which( ansDU$N.T < -0.2) ,]

ansDD = ans[ which( ans$BN.BT > 0.2) ,]
ansDD = ansDD[ which( ansDD$N.T > 0.2) ,]
c("up-up" , "upblood-downtissue" , "downblood-uptissue" , "down-down") -> s1
c(dim(ansUU)[1] , dim(ansUD)[1] , dim(ansDU)[1] , dim(ansDD)[1]) -> s2

