

address = "~/Desktop/newWork/R programming/data/homeWork1/fc7eeaab-cfbd-4e33-a8ab-154616aa8f68/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/TCGA-3N-A9WC-06A-11R-A38N-13.mirna.quantification.txt"
add = "~/Desktop/newWork/R programming/data/homeWork1/fc7eeaab-cfbd-4e33-a8ab-154616aa8f68/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3"


files = list.files(add)
files = files [-grep("isoform",files)]



sapply (files , FUN = TumororNormal) -> TorN_files

table_final <- lapply(c("T","N"), function(x) MakeTableFromFiles(files[TorN_files == x] , x))
do.call(cbind ,table_final ) -> table_final
write.table(table_final , add_wr , sep="\t")

data = read.table(address,sep="\t")


MakeTableFromFiles <- function(files , name ) {  
  #making table for tumor samples
  if (!length(files)) print ("Empty samples.")
  else {
    read =  lapply(paste ( add , files ,  sep = "/") , FUN = read.delim , header = TRUE) 
    do.call(cbind ,read ) -> table
    table <- table [,-grep("read_count|miRNA_ID|cross.mapped",colnames(table))]
    table <- log2 (table +1 )
    rownames (table )<- read [[1]]$miRNA_ID 
    setNames (table , paste0(name,1:(ncol(table))))
  }
}






TumororNormal <- function (str)
{
  check <- as.integer(substr(str,14,15))
  if ( check <10) return("T")
  else if (check >=10 && check <20 )return("N")
  else return("C")
}



#------------------------------------------------------
cancer_sample = sample(1:450,40)
sample_all_cancers = as.data.frame(table_final[[1]][,cancer_sample])


names(sample_all_cancers)[1] <- "x" 
do.call(vioplot,sample_all_cancers)

#--------------------------------------

cancer_samples = normalizeQuantiles(table_final[[1]])
normal_samples = normalizeQuantiles(table_final[[2]])

cancer_sample = sample(1:450,40)
sample_all_cancers = as.data.frame(cancer_samples[,cancer_sample])

names(sample_all_cancers)[1] <- "x" 
do.call(vioplot,c(sample_all_cancers,list(col="gold")))

final_data = cbind(cancer_samples,normal_samples)
ex = final_data
gr = substr(colnames(ex),1,1)

con = c("N-T")

answer = dif.ex(ex,gr,con)
#answer = answer[,c(2:15)]
answer = answer[order(answer$logFC,decreasing= T) , ]
