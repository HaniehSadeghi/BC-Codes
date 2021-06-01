setwd("~/Documents/Academy/Projects/TCGA/src")
source("~/Documents/Academy/R Scripts/AliLib.R")
source("~/Documents/Academy/R Scripts/Classify.R")
ai.algs <- c(svm=function(...)svm(kernel="linear",...))
ai.algs.pred="none"
options(mc.cores=detectCores())
logFC.threshold <- 2
dataFolder <- "../data/"
library(rgl)

MultiTissueClassify <- function(tissues, file.name="MultiTissue.pdf", only.upregulated.mirs = T) {
  x <- SelectFeatures(subset(mirs, tissue %in% tissues), "type", mirnames, 10, 10, 1)
  y <- lapply(x[1:5], function(a) setNames(a, gsub("hsa.mir.", "", names(a))))
  y <- lapply(y, function(a) data.frame(miRNAs = factor(names(a), levels=rev(names(a))), Accuracy = as.numeric(a)))
  y <- do.call(rbind, lapply(y, function(a)head(a,1)))
  y$miRNAs <- unlist(sapply(y$miRNAs, function(x) paste(collapse=", ", sort(unlist(strsplit(as.character(x), " "))))))
  y$miRNAs <- factor(y$miRNAs, levels=y$miRNAs)
  pdf(file.name,width=18)
  ggplot(y, aes(miRNAs, Accuracy,fill=miRNAs))+geom_bar(stat="identity")+theme_complete_bw() + theme(legend.position="none")+xlab("")+coord_flip()+scale_y_continuous(labels = percent)+geom_text(aes(label=round(Accuracy*100,0)))
  do()
  x
}

TissueClassify <- function(tissue, only.upregulated.mirs = T, iterations=1, ...) {
  mirs$Type <- ifelse(mirs$tissue==tissue, as.character(mirs$type), "Other")
  means <- aggregate(.~Type, mirs[,c(mirnames, "Type")], mean)
  rownames(means) <- means$Type
  x <- apply(mirs[,mirnames], 2, function(x) max(x) - min(x))
  mirnames <- mirnames[x>logFC.threshold]
  if (only.upregulated.mirs) mirnames <- mirnames[means["Tumor", mirnames] > means["Normal", mirnames] + 1]
  x <- SelectFeatures(mirs, "Type", mirnames, 5, 10, iterations, ...)
  y <- lapply(x[1:5], function(a) setNames(a, gsub("hsa.", "", gsub("mir.","", names(a), ignore.case = T))))
  y <- lapply(y, function(a) data.frame(miRNAs = factor(names(a), levels=rev(names(a))), Accuracy = as.numeric(a)))
  y <- do.call(rbind, lapply(y, function(a)head(a,1)))
  y$miRNAs <- unlist(sapply(y$miRNAs, function(x) paste(collapse=" ", sort(unlist(strsplit(as.character(x), " "))))))
  y$miRNAs <- factor(y$miRNAs, levels=y$miRNAs)
  pdf(paste0(tissue,".pdf"),width=18)
  print(ggplot(y, aes(miRNAs, Accuracy,fill=miRNAs))+geom_bar(stat="identity")+theme_complete_bw() + theme(legend.position="none",plot.background = element_rect(fill = '#FFE4E1'))+xlab("")+coord_flip()+scale_y_continuous(labels = percent)+geom_text(aes(label=round(Accuracy*100,0))))
  do()
  x
}

lapply(unique(mirs$tissue), TissueClassify, only.upregulated.mirs=F, iterations=1)

ReadData <- function(dataFolder) {
  ls <- list.files(dataFolder,"*.txt")
  mirs <- mclapply(paste0(dataFolder,ls), read.delim)
  mirs <- mclapply(mirs, function(x) x[sort(rownames(x)),])
  tissue <- rep(sub(".txt","",ls), times=sapply(mirs, ncol))
  type <- ifelse(substring(unlist(mclapply(mirs, colnames)),1,1)=="N","Normal","Tumor")
  mirs <- do.call(cbind, mirs)
  colnames(mirs)<-1:ncol(mirs)
  mirs <- data.frame(t(mirs))
  mirnas <- colnames(mirs)
  mirs$tissue <- factor(tissue)
  mirs$type <- factor(type)
  mirs$tt <- factor(paste(mirs$tissue, mirs$type,sep="_"))
  mirs <<- mirs
  mirnames <<- mirnas
}

ReadData(dataFolder)
mirs[,mirnames] <- t(normalizeQuantiles(t(mirs[,mirnames])))


Show3d <- function(tissues, x, y, z, radius=0.05, axesTitle=F) {
  type <- ifelse(mirs$tissue %in% tissues, "This", "Other")
  type <- paste0(type, mirs$type)
  col <- as.character(sapply(type, function(x) switch(x, ThisNormal="green", ThisTumor="red",OtherNormal="blue",OtherTumor="blue")))
  x <- paste0("hsa.miR.",x)
  y <- paste0("hsa.miR.",y)
  z <- paste0("hsa.miR.",z)
  open3d()
  spheres3d(mirs[,x], mirs[,y], mirs[,z], color=col, radius = ifelse(mirs$tissue%in%tissues, radius, radius/3))
  axes3d()
  if (axesTitle) title3d(xlab=x, ylab=y, zlab=z)
}

FixAnnotations <- function() {
  x <- read.delim("../ann/mimat_name.txt")
  
  y <- read.delim(paste0("../data/BLCA.txt"))
  y <- y[order(rownames(y)),]
  b <- rownames(y)
  a <- data.frame(full=b, ID=sub(".*\\|","",b), name=sub("\\|.*","",b))
  x <- x[,2:3]
  x <- x[!duplicated(x),]
  
  z <- merge(a, x, by.x="ID", by.y="Alias", all.x=T)
  z$Name <- as.character(z$Name)
  z[is.na(z$Name),"Name"] <- as.character(z[is.na(z$Name),"name"])
  rownames(z) <- z$full
  orn <- rownames(y)
  o <- order(-rowMeans(y))
  
  FixTissue <- function(tissue) {
    y <- read.delim(paste0("../data/", tissue,".txt"))
    y <- y[order(rownames(y)),]
    if (any(rownames(y)!=orn)) {
      cat("Error")
      return()
    }
    y <- y[o,]
    ynames <- z[rownames(y),"Name"]
    y <- y[!duplicated(ynames),]
    rownames(y) <- ynames[!duplicated(ynames)]
    write.table(y, file=paste0("../data/",tissue,"_2.txt"),row.names=T, col.names=T, sep="\t",quote=F)
  }
  
  invisible(lapply(unique(mirs$tissue), FixTissue))
}

