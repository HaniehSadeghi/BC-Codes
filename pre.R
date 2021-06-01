
library(pheatmap)
library(ggplot2)
library(reshape2)
library(grid)

add = "~/Desktop/Projects/Early_Detection_of_Cancer/Saarland_works/data/BRCA.txt"
data = read.delim(add)
data = normalizeQuantiles(data)
ma_plot <- function( data , name  , LogFC = 1.5 , adjpval = 0.001 , x = 1) { 
 

  ex = data
  gr = substr(colnames(ex),x,x)
  con = c("N-T")
  answer2 = dif.ex(ex,gr,con)
  ans2 = answer2[order(answer2$logFC,decreasing= F) , ]
  ans2 = ans2[ which ( ans2$adj.P.Val < adjpval) , ]
  ans2 = ans2[ which ( ans2$logFC > 1.5) , ]
  

  
  pdf(paste0( "~/Desktop/Saarland/results/MA/" , paste0(name , ".pdf") ),width=17,height=10)
  print(ggplot(ans2 , aes( x = AveExpr , y = logFC )) + 
    geom_abline(intercept  = 0 , slop = 0 , color = "red" , size = 2 , alpha = 3/4) +
    geom_abline(intercept  = LogFC , slop = 0 , color = "blue" , size = 1 , alpha = 1/2) +
    geom_abline(intercept  = -LogFC , slop = 0 , color = "blue" , size = 1 , alpha = 1/2) +
    geom_point(stat = "identity",size = 2) +
    theme(panel.background = element_rect(fill = 'white', colour = 'red'),
          panel.border = element_rect(fill = NA, colour = "black", size=2),
          strip.text.x = element_blank(),
          plot.background = element_rect(fill = '#FFE4E1'),
          strip.background = element_rect(fill = "#FFE4E1"),
          axis.text.x =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", vjust = 1),
          axis.text.y =       element_text(size = 40 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
          axis.title.x =  element_text(size = 40 * 0.8, lineheight = 0.9, colour = "black", vjust = 1),
          axis.title.y = element_text(size = 40 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
          legend.position = "none"
    ) + xlab("Average Expression") +
    ylab("Log Fold Change")  +
    ggtitle(name))
  dev.off()
}

rownames(ans)


old_data = read.delim("~/Desktop/newWork/data/ex_samples_tissues/BRCA.txt" , header = T)
old_data["hsa-mir-21",]
old_data_normal = old_data[,grep("N", colnames(old_data))]
mean(as.numeric(old_data_normal["hsa-mir-21",]))

data[grep("hsa-mir-21",rownames(data)),]
data_normal = data[,grep("N" , colnames(data))]
mean(as.numeric(data_normal["hsa-miR-21-5p",]))
mean(as.numeric(data_normal["hsa-miR-21-3p",]))


test = read.delim("~/Desktop/Saarland/TCGA data/BRCA/f7826293-2eb5-4928-9965-ddd81e76797e//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3//TCGA-E9-A1RH-11A-34R-A168-13.isoform.quantification.txt" , header = T)

test[which(test$miRNA_ID == "hsa-mir-21"),]
log2(sum(test[which(test$miRNA_ID == "hsa-mir-21"),4]))

tt = read.delim("~/Desktop/Saarland/TCGA data/BRCA/f7826293-2eb5-4928-9965-ddd81e76797e//miRNASeq//BCGSC__IlluminaHiSeq_miRNASeq//Level_3/TCGA-E9-A1RH-11A-34R-A168-13.mirna.quantification.txt" , header = T )
log2(tt[which(tt$miRNA_ID == "hsa-mir-21") ,3 ])




#################### 2D scatterplot BRCA 2 mirs
data[c("hsa-miR-21-5p","hsa-miR-183-5p"),]-> nn
nn = t(nn)
colnames(nn) <- c("R1", "R2")
nn = as.data.frame(nn)
Sample = substr(row.names(nn) , 1,1)
pdf(paste0( "~/Desktop/newWork/update/" , paste0("2D_scatterplot" , ".pdf") ),width=17,height=10)
print(ggplot(data = nn, aes(x = R1, y = R2, color=Sample)) +
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
              
        ) + xlab("hsa-miR-183-5p") +
          ylab("hsa-mir-21-5p") )

dev.off()

