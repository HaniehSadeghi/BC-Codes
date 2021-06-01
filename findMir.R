
# tis_names : esme hame bimari haye ke darim
# normals : sample haye normal tamami baft ha dar ghaleb yek table kenare yekdigar gharar begirand ( row = mirs , col = samples )
# library(limma)
# A = BLCA , B = BRCA , C = ESCA , D = HNSC , E = KICH , F = KIRC , G = KIRP , H = LIHC , I = LUAD , J = LUSC , K = PRAD , L = STAD , M = THCA , O = UCEC  
tis_names = c("BLCA" , "BRCA" , "ESCA" , "HNSC" , "KICH" , "KIRC" , "KIRP" , "LIHC" , "LUAD" , "LUSC" , "PRAD" , "STAD" , "THCA" , "UCEC" )
normals = NormalsTable();
cancers = CancerTable();
tt = cbind(cancers,normals)
tt = normalizeQuantiles(tt)
normals = tt[,grep('N' , colnames(tt))]
cancers = tt[,grep('T' , colnames(tt))]


# in dictionary niaz mishe toye kar
myMap = vector(mode="list" , length= 14)
names(myMap) <- tis_names
myMap[[1]] = "A";myMap[[2]] = "B";myMap[[3]] = "C";myMap[[4]] = "D";myMap[[5]] = "E";myMap[[6]] = "F";myMap[[7]] = "G";
myMap[[8]] = "H";myMap[[9]] = "I";myMap[[10]] = "J";myMap[[11]] = "K";myMap[[12]] = "L";myMap[[13]] = "M";myMap[[14]] = "O";

# joda joda ba har baft baresi kone ejtema begire 
# name : esme saratan X + .txt
# negative = TRUE bashe mir haye O mide ke dar saratani bayaneshun bisthare az kolle normal ha ,
# negative = FALSE bashe mir haye O mide ke dar saratani bayane kamtari daran nesbat be kolle normal
# num = tedad mir haye ke bad az har dif.ex intersect migire 

# bayad faile haye ex_sample_tissues load shode bashad






# in ye jadval misaze az hameye normal ha va saratani baft X va dif ex migire Y taye avalO barmigardune 
# tis_names  = name hameye bimari ha 
# num = tedad mir haye ke dar akhar baramun dar khoroji mide
# name = X yek baft khas ke bayad 
# LogFC = -1.5 . anhaye ra peyda mikone ke logFC kamtar az in bashe
#  adjpval = 0.001 anhaye ra peyda mikone ke in shart ro hamdashte bashan


topMirs_AllNomrals<- function( name , LogFC = -1.5 , adjpval = 0.001 ) {
  
  

  # hazf bimari name az list 
  index = which(tis_names == name)
  tmp = tis_names[-index]
  
  # bimari mored nazar Ro mikhune 
  tmp = colnames(cancers)
  grep(myMap[[name]] , tmp) -> indexes
  cancer = cancers[,indexes]
  colnames(cancer)<- gsub( myMap[[name]] , "", colnames(cancer)) 
  
  tmp2 = normals
  colnames(tmp2) <- sapply(1:dim(normals)[2] , function(x) { paste("N" , as.character(x) , sep ="")} )
  data = cbind(tmp2,cancer)
  ex = data
  gr = substr(colnames(data),1,1)
  con = c("N-T")
  answer = dif.ex(ex,gr,con)
  answer = answer[order(answer$logFC,decreasing= F) , ]
  ans = answer[ which ( answer$logFC < LogFC) , ]
  ans = ans[ which ( ans$adj.P.Val < adjpval) , ]
  rownames(ans)
}



# general 
# ye table dorost mikone az hame normal ha 
# updated with normalization
NormalsTable <-function (){
    
  add = "~/Desktop/newWork/data/ex_samples_tissues/"
  files = list.files(add)
  tables = lapply( paste(add , files , sep="/" ) , FUN=read.delim , header = T ) 
  cl = c( "A" , "B" , "C" , "D" , "E" , "F" , "G" , "H" , "I" ,"J" , "K" , "L" , "M" , "O" )
  for ( i in 1:14 ) { colnames(tables[[i]]) = paste( cl[i] , colnames(tables[[i]]) , sep ="" ) } 
  do.call ( cbind , tables ) -> all_mirs_exp  
  colnames(all_mirs_exp) -> cnames  
  normal_indexes = grep('N' , cnames)
  all_normals = all_mirs_exp[,normal_indexes]
  normalized_normals = normalizeQuantiles(all_normals)
  normalized_normals
}


CancerTable <-function (){
  
  add = "~/Desktop/newWork/data/ex_samples_tissues/"
  files = list.files(add)
  tables = lapply( paste(add , files , sep="/" ) , FUN=read.delim , header = T ) 
  cl = c( "A" , "B" , "C" , "D" , "E" , "F" , "G" , "H" , "I" ,"J" , "K" , "L" , "M" , "O" )
  for ( i in 1:14 ) { colnames(tables[[i]]) = paste( cl[i] , colnames(tables[[i]]) , sep ="" ) } 
  do.call ( cbind , tables ) -> all_mirs_exp  
  colnames(all_mirs_exp) -> cnames  
  cancer_indexes = grep('T' , cnames)
  all_cancers = all_mirs_exp[,cancer_indexes]
  normalized_cancers = normalizeQuantiles(all_cancers)
  normalized_cancers

}

# just write files on mirs folder
write_names <- function( list , name ) {
  add = "/Users/Aryan/Desktop/newWork/Mirs/";
  add = paste ( add , name , sep = "");
  
  write.table(list , add , sep = "\t")
}






# need to be update 

topMirs_Split<- function( name , num = 50 , negative = TRUE) { 
  
  address = "~/Desktop/newWork/data/ex_samples_tissues/"
  address = paste(address , name , sep = "");
  table = read.table(address , sep = "\t" );
  tmp = colnames(table)
  grep('T' , tmp) -> indexes
  cancer = table[,indexes]
  ans = rownames(cancer)
  
  for ( i in 1:14 ) {
    address = "~/Desktop/newWork/data/ex_samples_tissues/"
    address = paste(address , paste ( tis_names[i] , ".txt" , sep ="") , sep ="");
    tab = read.table(address , sep = "\t" );
    tmp = colnames(tab)
    grep('N' , tmp) -> indexes
    norm = tab[,indexes]
    
    data = cbind( norm , cancer ) 
    ex = data
    gr = substr(colnames(data),1,1)
    con = c("N-T")
    answer = dif.ex(ex,gr,con)
    if ( negative ) {
      answer = answer[order(answer$logFC,decreasing= F) , ]
      ans2 = rownames( head(answer,num))
      ans = intersect(ans , ans2)
    } else {
      answer = answer[order(answer$logFC,decreasing= T) , ]
      ans2 = rownames( head(answer,num))
      ans = intersect(ans , ans2)
    }
    
  }
  return(ans)
  
}