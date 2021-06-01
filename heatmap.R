# heat mape va hala nemudar haye dge mikeshim abraye mir hayee ke cadidate shodan baraye har bimari 


#heatmap mikeshe baraye har bimari ke tasire mir haye X toye in bimari nesbat be baghiye chejuriyan

# name = BRCA , esme bimari ke default breast cancer hast 
# num = tedad mir ha
# data = ye file ke satrash mir ha hast + soton hash bimari ha var har khune dif.ex hast
# sourceAddress = jayee ke nemudar ha keshide mishan va unja gharar migiran
# FileName = esme fileE ke save mikone 


# ina niaz mishe to tabe mirs_data 
normals = NormalsTable();
cancers = CancerTable();
myMap = vector(mode="list" , length= 14)
names(myMap) <- tis_names
myMap[[1]] = "A";myMap[[2]] = "B";myMap[[3]] = "C";myMap[[4]] = "D";myMap[[5]] = "E";myMap[[6]] = "F";myMap[[7]] = "G";
myMap[[8]] = "H";myMap[[9]] = "I";myMap[[10]] = "J";myMap[[11]] = "K";myMap[[12]] = "L";myMap[[13]] = "M";myMap[[14]] = "O";



HeatMap <- function( name = "BRCA" , SourceAddress = "/Users/Aryan/Desktop/newWork/graphs/" , FileName = "heatMap.pdf" ){
    list = topMirs_AllNomrals(name=name)
 
    data = mirs_data(name)
    my_data = data[list,]
    setwd(SourceAddress)
    pdf(FileName, width=14, height=14)
    green = colorRampPalette(c("green","black"))(128)
    red = colorRampPalette(c("black","red"))(128)
    bk = unique(c(seq(-10, 0, length=127), seq(0,7 , length=128)))
    pheatmap(my_data, border_color=NA, color=c(green , red) , breaks=bk)
    dev.off()
}

# nemitone vaghti run mikonam khub bekeshe dasti bayad bekesham
PCA <- function( name = "BRCA" , SourceAddress = "/Users/Aryan/Desktop/newWork/graphs/" , FileName = "heatMap.pdf" ) {
  
  list = topMirs_AllNomrals(name=name)
  
  data = mirs_data(name)
  my_data = data[list,]
  
  pc = prcomp(my_data)
  scores = as.data.frame(pc$x)
  
  
 
  
  setwd(SourceAddress)
  pdf(FileName, width=20, height=14)

  ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_text(colour = "tomato", alpha = 0.8, size = 4) +
    ggtitle(  paste( "PCA plot for" , name , sep = " ") )  
  
  
  dev.off()
}




# in tabe miad baraye bimari X cancer O var midare varoye tak take bimari haye dge for mizane va dif.ex hessab mikone 
# va dar ghaleb ye table mide be ma 
# satr ha = mir ha 
# soton ha = bimari 
# a(i,j) = mir i am , in mir age sarataniye az X ro dashte bashim va normal az J taghiratesh chejuriye

mirs_data <- function( name = "BRCA", negative = TRUE) {
  
  tmp = colnames(cancers)
  grep(myMap[[name]] , tmp) -> indexes
  cancer = cancers[,indexes]
  colnames(cancer)<- gsub( myMap[[name]] , "", colnames(cancer)) 
  
  my_ans = c()
  my_ans = cbind(my_ans,rownames(cancer))
  rownames(my_ans) <- rownames(cancer)
  col = c("V1")
  
  
  for ( i in 1:14 ) {
    
    tmpName = tis_names[i]
    tmp = colnames(normals)
    grep(myMap[[tmpName]] , tmp) -> indexes
    norm = normals[,indexes]
    colnames(norm)<- gsub( myMap[[tmpName]] , "", colnames(norm)) 
    
    data = cbind( norm , cancer ) 
    ex = data
    gr = substr(colnames(data),1,1)
    con = c("N-T")
    answer = dif.ex(ex,gr,con)
    
    
    answer2 = cbind(rownames(answer) , answer)
    colnames(answer2) -> ll
    ll[1] = "V1"
    colnames(answer2) <- ll

    my_ans = merge( my_ans, answer2 , by = "V1")
    tmp = i + 1 
    my_ans = my_ans[,c(1:tmp)]
    col = c(col , tis_names[i])
    colnames(my_ans) = col
    rownames(my_ans) = my_ans[,"V1"]
      
    
  }
  my_ans = my_ans[,-1]
  return(my_ans)
}