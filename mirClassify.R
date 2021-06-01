#source("~/Documents/Academy/R Scripts/AliLib.R")
#source("~/Documents/Academy/R Scripts/Classify.R")
ai.algs <- c(svm=svm)
ai.algs.pred="none"
options(mc.cores=detectCores())

logFC.threshold <- 2
dataFolder <- "../Early_Detection_of_Cancer/Saarland_works/data/"
ls <- list.files(dataFolder,"*.txt")
mirs <- mclapply(paste0(dataFolder,ls), read.delim)
mirs <- mclapply(mirs, function(x) x[sort(rownames(x)),])
tissue <- rep(sub(".txt","",ls), times=sapply(mirs, ncol))
type <- substring(unlist(mclapply(mirs, colnames)),1,1)
mirs <- do.call(cbind, mirs)
colnames(mirs)<-1:ncol(mirs)
mirs <- data.frame(t(mirs))
mirnas <- colnames(mirs)[apply(mirs, 2, max) - apply(mirs, 2, min) > logFC.threshold]
mirs$tissue <- factor(tissue)
mirs$type <- factor(type)
mirs$tt <- factor(paste(mirs$tissue, mirs$type,sep="_"))
Classify(mirs, "tt", mirnas, alg = svm, name = "svm",pred.type = "none") -> tmp_ml



##### just BRCA 
mirs_BRCA = mirs[mirs$tissue == "BRCA",]
comb1_BRCA = ClassifyCombinations(mirs_BRCA,alg = svm, "type", mirnas, comb.num = 1, iterations=10)
comb2_BRCA = ClassifyCombinations(mirs_BRCA, "type", head(names(comb1_BRCA),20), 2, 10)
comb3_BRCA = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb2_BRCA)," "))), 3, 10)
comb4_BRCA = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(head(names(comb3_BRCA),10)," "))), 4, 10)
comb5_BRCA = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb4_BRCA)," "))), 5, 10)
comb6_BRCA = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb5_BRCA)," "))), 6, 10)

#### try with cross validation 
mirs_BRCA = mirs[mirs$tissue == "BRCA",]
comb1_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", mirnas, 1, 4)
comb2_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", head(names(comb1_BRCA_cr),20), 2, 10)
comb3_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb2_BRCA_cr)," "))), 3, 10)
comb4_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(head(names(comb3_BRCA_cr),10)," "))), 4, 10)
comb5_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb4_BRCA_cr)," "))), 5, 10)
comb6_BRCA_cr = ClassifyCombinations(mirs_BRCA, "type", unique(unlist(strsplit(names(comb5_BRCA_cr)," "))), 6, 10)



#### I have added this part for 3 groups : BRCA_T Normal and  other cancers 
mirs$type_3group = rep("0" , nrow(mirs))
unique(as.character(mirs$tt[grep("_N" , mirs$tt)])) -> all_normal
unique(as.character(mirs$tt[grep("_T" , mirs$tt)])) -> all_cancer
all_cancer_BRCA = all_cancer[ all_cancer != "BRCA_T"]
BRCA = "BRCA_T"

mirs$type_3group[mirs$tt %in% all_normal ] = "N"
mirs$type_3group[mirs$tt %in% all_cancer_BRCA ] = "T"
mirs$type_3group[mirs$tt %in% BRCA ] = "BR"

####### 
comb1 = ClassifyCombinations(mirs, "type_3group", mirnas, 1, 4)
comb2 = ClassifyCombinations(mirs, "type_3group", head(names(comb1),20), 2, 10)
comb3 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb2)," "))), 3, 10)
comb4 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(head(names(comb3),10)," "))), 4, 10)
comb5 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb4)," "))), 5, 10)
comb6 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb5)," "))), 6, 10)

comb7 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb6)," "))), 7, 10)
comb8 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb7)," "))), 8, 10)
comb9 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb8)," "))), 9, 10)
comb10 = ClassifyCombinations(mirs, "type_3group", unique(unlist(strsplit(names(comb9)," "))), 10, 10)


mirs$Tumor_Type <- as.character(mirs$tissue)
mirs$Tumor_Type[mirs$type=="Normal"]= "Normal"

pdf("~/Desktop/CancerMir.pdf",width=10,height=10)
ggplot(mirs, aes(hsa.mir.183, hsa.mir.21, color=Tumor_Type))+geom_point(size=1.5)+geom_dl(aes(label=Tumor_Type, color=Tumor_Type),method="smart.grid")+theme_complete_bw()+theme(legend.position="none")
dev.off()

pdf("~/Desktop/CancerMirTissue.pdf",width=30,height=30)
ggplot(mirs, aes(hsa.mir.183, hsa.mir.21, color=type))+geom_point(size=3)+theme_complete_bw()+facet_wrap(~tissue)
ggplot(mirs, aes(hsa.mir.183, hsa.mir.21, color=type))+geom_point(size=3)+theme_complete_bw()+facet_wrap(~tissue,scales = "free")
dev.off()
ggplot(mirs, aes(hsa.let.7a.1, hsa.let.7b, color=type))+geom_point(size=3)+theme_complete_bw()+facet_wrap(~tissue,scales = "free")

mirs$Brca <- ifelse(mirs$tissue=="BRCA",mirs$type,"Other")
ggplot(mirs, aes(hsa.mir.10b, hsa.mir.190b, color=Brca))+geom_point(size=4)+theme_complete_bw()

mirs<- mirs[,c("hsa.mir.139", "hsa.mir.199a.2", "hsa.mir.27b")]



##### new task 
sub_mirs = c("hsa.miR.10b.5p" , "hsa.miR.193a.5p" , "hsa.miR.190b" , "hsa.miR.196a.5p" , "hsa.miR.199b.5p")
sub_mirs2 = c("hsa.miR.10b.3p" ,  "hsa.miR.193a.5p"  , "hsa.miR.190b" , "hsa.miR.196a.5p" , "hsa.miR.21.5p")
res = Classify(mirs,target.col = "type_3group",trainSize=19*nrow(mirs)/20,testSize=nrow(bm)/20,cols = sub_mirs ,quiet = T)
res
f = -1
for ( i in 1:100) {
  res = Classify(mirs,target.col = "type_3group",trainSize=19*nrow(mirs)/20,testSize=nrow(bm)/20,cols = sub_mirs ,quiet = T)
  f = max (f, res$accuracy)
}
f
