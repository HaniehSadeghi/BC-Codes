library (limma)
library (parallel)
addBRCA =  "~/Desktop/Saarland/TCGA data/BRCA/f7826293-2eb5-4928-9965-ddd81e76797e//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addBLCA = "~/Desktop/Saarland/TCGA data/BLCA/bc0f8364-55ed-4778-b7d8-6443e44613c8//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addTHCA = "~/Desktop/Saarland/TCGA data/THCA/57ccbe9f-586c-452c-90f1-9b045619008a//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addESCA = "~/Desktop/Saarland/TCGA data/ESCA/01d0c1c7-26a8-4615-a7a0-4ca7d7811d3b//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addHNSC = "~/Desktop/Saarland/TCGA data/HNSC/29504953-d519-4429-be37-71436397786d//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addKICH = "~/Desktop/Saarland/TCGA data/KICH/bd357678-e88e-40a2-9777-7957d46902c3//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addKIRC = "~/Desktop/Saarland/TCGA data/KIRC/24df071a-697b-401e-833a-4bed018fb94d//miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addKIRP = "~/Desktop/Saarland/TCGA data/KIRP/158af0f4-8da4-42cd-98c4-f6e63ed55acd//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addLIHC = "~/Desktop/Saarland/TCGA data/LIHC/c10ec660-8ede-485d-a018-93ece0695313//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addLUAD = "~/Desktop/Saarland/TCGA data/LUAD/2c90338f-b8be-40e0-8c00-2ca9d34b653c//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addPRAD = "~/Desktop/Saarland/TCGA data/PRAD/671c54a5-bff2-4d0c-b24d-cf1968609d7c//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addSTAD = "~/Desktop/Saarland/TCGA data/STAD/9ade7af0-e141-44ed-9ded-4603120b07ee//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addUCEC = "~/Desktop/Saarland/TCGA data/UCEC/2535c385-cbd1-4ab5-8dc3-a2f9e5c67331//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
addLUSC = "~/Desktop/Saarland/TCGA data/LUSC/bda547a7-fc86-4521-9822-b48f72bc919a//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/"
add_wr  = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/UCEC.txt"
addBRCAwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/BRCA.txt"
addBLCAwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/BLCA.txt"
addTHCAwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/THCA.txt"
addESCAwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/ESCA.txt"
addHNSCwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/HNSC.txt"
addKICHwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/KICH.txt"

addKIRCwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/KIRC.txt"
addKIRPwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/KIRP.txt"
addLIHCwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/LIHC.txt"
addLUADwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/LUAD.txt"
addPRADwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/PRAD.txt"
addSTADwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/STAD.txt"
addUCECwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/UCEC.txt" 
addLUSCwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/LUSC.txt"

add = c(addBRCA,addBLCA,addTHCA,addESCA,addHNSC,addKICH,addKIRC,addKIRP,addLIHC,addLUAD,addLUSC,addPRAD,addSTAD,addUCEC)
addwr = c(addBRCAwr,addBLCAwr,addTHCAwr,addESCAwr,addHNSCwr,addKICHwr,addKIRCwr,addKIRPwr,addLIHCwr,addLUADwr,addLUSCwr,addPRADwr,addSTADwr,addUCECwr)
name2 = c("brca" , "blca" , "thca" , "esca" , "hnsc" , "kich" , "kirc" , "kirp" , "lihc" ,"luad" ,"lusc" , "prad" ,"stad" , "ucec")

# Reads a group of files, return in a table the "reads_per_million_mapped columns of 
#    TCGA data

MakeTableFromFiles <- function(files , name ) {  
  #making table for tumor samples
  if (!length(files)) print ("Empty samples.")
  else {
    read =  mclapply(paste ( add , files ,  sep = "/") , FUN = read.delim , header = TRUE) 
    do.call(cbind ,read ) -> table
    table <- table [,-grep("read_count|miRNA_ID|cross.mapped",colnames(table))]
    table <- log2 (table +1 )
    rownames (table )<- read [[1]]$miRNA_ID
    setNames (table , paste0(name,1:(ncol(table))))
  }
}
# dat is for isoforms
MakeTableFromFiles2 <- function(add , files , name ) {  
  #making table for tumor samples
  if (!length(files)) print ("Empty samples.")
  else {
    rows = c()
    read =  mclapply(paste ( add , files ,  sep = "/") , FUN = read.delim , header = TRUE) 
    # calculate read counts 
    
    for ( i in 1:length(read) ) {
      read[[i]] = makeTableFromIso(read[[i]])
      rows = unique(c(rows,read[[i]]$miRNA_ID))
    }
   
    for ( i in 1:length(read) ) {
      read[[i]] = addZero ( read[[i]] , rows)
      
    }
    do.call(cbind ,read ) -> table
    table <- table [,-grep("miRNA_ID",colnames(table))]
    rownames (table)<- read [[1]]$miRNA_ID
   
    setNames (table , paste0(name,1:(ncol(table))))
  }
}

makeTableFromIso <- function(table) {

  table = table[,c(1,4,6)]
  table = table[-c ( which(table$miRNA_region == "precursor") , which(table$miRNA_region == "unannotated") , which(table$miRNA_region == "stemloop")),]  
  table$miRNA_region = sapply(1:length(table$miRNA_region) , function(x) { strsplit2(table$miRNA_region[x], split=",")[2]})
  table$miRNA_ID = paste( table$miRNA_ID , table$miRNA_region , sep = "|" ) 
  table$miRNA_region = NULL
  table = aggregate( . ~ miRNA_ID, data = table, sum) 
  table$exp = log2(table$reads_per_million_miRNA_mapped + 1 )
  table$reads_per_million_miRNA_mapped = NULL
  table
}

addZero <- function(table, rows) {
  df = data.frame(miRNA_ID = rows, exp = rep(0,length(rows)))
  
  rbind(df , table) -> table
  table = aggregate( . ~ miRNA_ID, data = table, sum) 
  table
  
}

addZero2 <- function(table, rows) { 
  adder = rep(0,dim(table)[2])
  for ( i in 1:length(rows) ) {
    if(is.na(table[rows[i],1])){
      table = rbind(table , adder)
      rownames(table)[dim(table)[1]] <- rows[i] 
    }
  }
  table
}


readCount <- function(add , name){
  files = list.files(add)
  iso_files = files [ grep("isoform", files)]
  sapply (iso_files , FUN = TumororNormal) -> TorN_files
  files = iso_files
  read =  mclapply(paste ( add , files ,  sep = "/") , FUN = read.delim , header = TRUE) 
  rcounts = sapply(1:length(read) , function(x) { sum(read[[x]]$read_count)})
  which((TorN_files == "N") == TRUE ) -> nor
  which((TorN_files == "N") == FALSE ) -> tur
  rcounts[tur] ->rct
  rcounts[nor] ->rcn
  tumor = c(mean(rct) , var(rct)  , sd(rct) )
  normal = c(mean(rcn) , var(rcn) , sd(rcn) )
  df = t(data.frame(tumor,normal))
  colnames(df) <- c("Average" , "Variance" , "Standard Deviation")
  rownames(df) <- c(paste("tumor" , name , sep="_") ,paste("normal" , name , sep="_") )
  df
}
#main function,
#this function make table for samples in a tissue,
#columns are miRNA expression in samples, rownames are miRNA_ID. 
#normal samples in one group, tumor samples in one other group
#add = the address of samples
#add_wr = the address to save the table as the output


#add="/home/a.sharifi/TCGAmir/Data/BRCA/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3"
#add_wr = "/home/a.sharifi/TCGAmir/Data/compiled_data/diff_ex_tissues"

# returns #of normal samples and #of tumor samples
# table E ke write mikone b onvane khoruji, mishe vorudie dif.ex
makeTableforTCGAsamples_forTissue <- function (add3, add_wr)
{
  files = list.files(add3)
  iso_files = files [ grep("isoform", files)]
  sapply (iso_files , FUN = TumororNormal) -> TorN_files
  
  table_final <- lapply(c("T","N"), function(x) MakeTableFromFiles2(add3 , iso_files[TorN_files == x] , x))
  if ( FALSE ) {
    table_final[[1]]$name = rownames(table_final[[1]])
    table_final[[2]]$name = rownames(table_final[[2]])
    t = merge(table_final[[1]] , table_final[[2]] , by.x = "name" , by.y = "name")
    t$name = NULL
    myData = t;
  }
  
  rows = unique( c(row.names(table_final[[1]]) , rownames(table_final[[2]])))
  table_final1 = addZero2(table_final[[1]],rows)
  table_final2 = addZero2(table_final[[2]],rows) 
  table_final1$name = rownames(table_final1)
  table_final2$name = rownames(table_final2)
  table_final1 = table_final1[ order(row.names(table_final1)) , ]
  table_final2 = table_final2[ order(row.names(table_final2)) , ]
  t = merge(table_final1,table_final2 , by.x = "name" , by.y = "name" )
  row.names(t) <- t$name
  t$name = NULL
  write.table(t , add_wr , sep="\t")
  
  if(FALSE) {
    files = files [-grep("isoform",files)]
    #seperating normal samples from tumor samples
    sapply (files , FUN = TumororNormal) -> TorN_files
    #making table with normal samples in first columns and tumor samples in the last columns
    table_final <- lapply(c("T","N"), function(x) MakeTableFromFiles(files[TorN_files == x] , x))
    do.call(cbind ,table_final ) -> table_final
    write.table(table_final , add_wr , sep="\t")
  }
  #first is normal count
}  

#### Differential Expression
# ex: A matrix of gene expression values, genes in rows, samples in columns
# gr: A vector of factors/characters/integers of size equal to ncol(ex) indicating the cell type each sampleb belongs to
# con: The contrast of differential expression: c("gr1-gr2",...)
# ann: Annotation columns for the genes
# number: Maximum number of differentially expressed genes to report
dif.ex <- function(ex, gr, con, ann=NULL, number=Inf) {
  gr <- factor(gr)
  mat =data.frame(G=gr)
  rownames(mat) <- colnames(ex)
  design = model.matrix(~G+0,mat)
  colnames(design) <- levels(gr)
  fit <- lmFit(ex, design)
  cont.matrix <- makeContrasts(contrasts=con, levels=design)
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT=topTable(fit2, adjust="fdr", sort.by="B", number=number )
  if (! is.null(ann)) tT <- cbind(tT, ann[rownames(tT),])
  tT
}


#Tumor types range from 01 - 09, Normal types from 10 - 19, sample type=[14,15]
#this function specify the health of sample base on sample type.
TumororNormal <- function (str)
{
  check <- as.integer(substr(str,14,15))
  if ( check <10) return("T")
  else if (check >=10 && check <20 )return("N")
  else return("C")
}


Run_dif.ex <- function( add , writeAdd ) {
  ex = read.delim(add , header = T );
  gr = substr (colnames(ex) , 1 , 1)
  con = c("N-T");
  answer = dif.ex(ex,gr,con);
  write.table(answer , writeAdd , sep ="\t");
}


######### different sample miRNAs ro mikham yeki konam 
######### data haro download konam ye rows besazam bad baraye hamashun inO chap konam dobare 

add2 = "~/Desktop/Saarland/TCGA data/ex_samples_tissues/"
addwr = "~/Desktop/Saarland/TCGA data/ex_samples_tissues2/"
files = list.files(add2)
table =  lapply(paste ( add2, files ,  sep = "/") , FUN = read.delim , header = TRUE)
rows = c()
for ( i in 1:14 ) {
rows = unique( c(rows, rownames(table[[i]])) ) 
}

for ( i in 1:14 ) {
  table[[i]] = addZero2(table[[i]] , rows)
}

for ( i in 1:14 ) {
  write.table(table[[i]] , paste ( addwr, files[i] ,  sep = "/") , sep="\t")
}
