# graph for data sets 

library(reshape2)
library(dplyr)
library(ggplot2)
samples = c( 285+19 , 9+72 , 45+331 , 38+276 , 59+507 , 21+393 , 25+66 , 32+198 , 50+330 , 87+695, 71+237 , 46+436, 43+420,50+140)
cancer = c( 285 , 72 , 331 , 276 , 507 , 393 , 66 , 198 , 330 , 695 , 237 , 436 , 420 , 140)
normal = c(19 , 9 , 45 , 38 , 59 , 21 , 25 , 32 , 50 , 87 , 71 , 46 , 43 , 50)



############## barplot for sample size 

data = rbind(cancer ,normal)
data = rbind( data , samples)
data = as.data.frame(data)
colnames(data) <- c("BLCA" , "ESCA" , "LUSC" , "STAD" , "THCA" , "UCEC" ,"KICH" ,"KIRP" , "PRAD" , "BRCA" ,"KIRC" ,"LUAD" , "HNSC" ,"LICH" )

data = data[, order(data[3,],decreasing= T) ]
data = t(data)
data = as.data.frame(data)
data$name = rownames(data)
rownames(data) <- NULL

data <- melt(data,id.vars = "name")
data = data[1:28,]
colnames(data)<-c("Tissue_name" , "Type", "Number_of_Samples")
addGraph = "~/Desktop/newWork/new version of graphs/barplot_samplesSize.pdf"

pdf(addGraph,width=17,height=10)
ggplot(data , aes(x= Tissue_name , y = Number_of_Samples ,fill=Type )) + 
  geom_bar( stat="identity")+ 
  geom_bar(stat="identity",colour="black", show_guide=FALSE)+
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.border = element_rect(fill = NA, colour = "black", size=2),
        plot.background = element_rect(fill = "#FFE4E1"),
        strip.text.x = element_text(size = 25, colour = "black"),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(colour="white"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(1.5, "lines"),
        legend.text =       element_text(size = 20 * 1),
        legend.title =      element_text(size = 20 * 1),
        legend.position =   c(0.8,0.90)
      ) + xlab("Tissue Name") +
          ylab("Number of Sample") 
dev.off()


################## barplot for sample size END 

################### pie charts 

####### only for number of samples


Type = c("Cancer" , "Normal")
numSam = c(5493,611)
df <- data.frame(numSam , Type)
df <- df %>% mutate(pos = cumsum(numSam) -  numSam/2)
addGraph = "~/Desktop/newWork/new version of graphs/SampleSize.pdf"

pdf(addGraph,width=17,height=10)
ggplot(df , aes(x="", y=numSam,fill=Type )) + 
  geom_bar( width=1 ,stat="identity")+ 
  geom_text(aes( y=pos,label=numSam), size=25) +  # note y = pos
  coord_polar("y") +
  theme(panel.background = element_rect(fill = '#FFE4E1', colour = 'white'),
        panel.border = element_rect(fill = NA, colour = "white", size=0),
        strip.text.x = element_text(size = 25, colour = "black"),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(fill = "#FFE4E1",colour="white"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(5.5, "lines"),
        legend.text =       element_text(size = 20 * 2),
        legend.title =      element_text(size = 20 * 2),
        legend.position =   "right",
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) + xlab("") +
  ylab("") 
dev.off()
####### only for number of samples END

Type = c("BRCA vs All Normals" , "All Cancers vs All Normals")
numSam = c(54,27)
df2 <- data.frame(numSam , Type)
df2 <- df2 %>% mutate(pos = cumsum(numSam) -  numSam/2)

pdf("~/Desktop/newWork/slide/allcancers.pdf",width=17,height=10)
ggplot(df2 , aes(x="", y=numSam,fill=Type )) + 
  geom_bar( width=1 ,stat="identity")+ 
  geom_text(aes( y=pos,label=numSam), size=15) +  # note y = pos
  coord_polar("y") +
  scale_fill_manual(values=c("#F8766D","#00bfc4"))+
  theme(panel.background = element_rect(fill = '#FFE4E1', colour = 'white'),
        panel.border = element_blank(),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.text.x = element_text(size = 25, colour = "black"),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(fill ="#FFE4E1" , colour="#FFE4E1"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(3.5, "lines"),
        legend.text =       element_text(size = 40 ),
        legend.title =      element_text(size = 40 ),
        legend.position =   "right",
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) + xlab("") +
  ylab("") 
dev.off()



################## for all other feature dat we selecet some mirs 
Type = c("Less than 5" , "More than 5")
numSam = c(22,34)
df2 <- data.frame(numSam , Type)
df2 <- df2 %>% mutate(pos = cumsum(numSam) -  numSam/2)

pdf("~/Desktop/newWork/slide/mean.pdf",width=17,height=10)
ggplot(df2 , aes(x="", y=numSam,fill=Type )) + 
  geom_bar( width=1 ,stat="identity")+ 
  geom_text(aes( y=pos,label=numSam), size=15) +  # note y = pos
  coord_polar("y") +
  scale_fill_manual(values=c("#008B00","#FF6103"))+
  theme(panel.background = element_rect(fill = '#FFE4E1', colour = 'white'),
        panel.border = element_blank(),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.text.x = element_text(size = 25, colour = "black"),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(fill ="#FFE4E1" , colour="#FFE4E1"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(5.5, "lines"),
        legend.text =       element_text(size = 30 ),
        legend.title =      element_text(size = 30 ),
        legend.position =   "right",
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 45, colour = "black")
  ) + xlab("") +
  ylab("") + labs(title = "Mean")
dev.off()



#############################
Type = c("not Reported","Seretion Reported")
numSam = c(1,53)
df2 <- data.frame(numSam , Type)
df2 <- df2 %>% mutate(pos = cumsum(numSam) -  numSam/2)

pdf("~/Desktop/newWork/slide/secretion.pdf",width=17,height=10)
ggplot(df2 , aes(x="", y=numSam,fill=Type )) + 
  geom_bar( width=1 ,stat="identity")+ 
  geom_text(aes( y=pos,label=numSam), size=15) +  # note y = pos
  coord_polar("y") +
  scale_fill_manual(values=c("#F8766D","#00bfc4"))+
  theme(panel.background = element_rect(fill = '#FFE4E1', colour = 'white'),
        panel.border = element_blank(),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.text.x = element_text(size = 25, colour = "black"),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(fill ="#FFE4E1" , colour="#FFE4E1"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(5.5, "lines"),
        legend.text =       element_text(size = 30 ),
        legend.title =      element_text(size = 30 ),
        legend.position =   "right",
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 45, colour = "black")
  ) + xlab("") +
  ylab("") + labs(title = "Secretion")
dev.off()


############## chromosome
Type = c("2","5","12","15","9","3","17","19","11","14","16","1","7","X")
numSam = c(3,2,6,3,4,2,4,2,4,5,3,9,3,4)
df2 <- data.frame(numSam , Type)
df2 <- df2 %>% mutate(pos = cumsum(numSam) -  numSam/2)
col = c(7,10,12,28,21,32,33,31,53,84,83,101,124,154)

pdf("~/Desktop/newWork/slide/previous.pdf",width=17,height=10)
ggplot(df2 , aes(x="", y=numSam,fill=Type )) + 
  geom_bar( width=1 ,stat="identity")+ 
  geom_text(aes( y=pos,label=numSam), size=15) +  # note y = pos
  coord_polar("y") +  scale_fill_manual(values=color[col])+
  theme(panel.background = element_rect(fill = '#FFE4E1', colour = 'white'),
        panel.border = element_blank(),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.text.x = element_text(size = 25, colour = "black"),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.background = element_rect(fill ="#FFE4E1" , colour="#FFE4E1"), 
        legend.key =        element_rect(fill ="#FFE4E1", colour = "white", size = 0.25 ),
        legend.key.size =   unit(2.5, "lines"),
        legend.text =       element_text(size = 30 ),
        legend.title =      element_text(size = 30 ),
        legend.position =   "right",
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 45, colour = "black")
  ) + xlab("") +
  ylab("") + labs(title = "Chromosome")
dev.off()





################ barplot for logFC



# UP va DOWN bayad load beshan va hamchenin all2 ke mishe kolle mir ha

all = cbind(UP,DOWN)


UP$name = rownames(UP)
DOWN$name = rownames(DOWN)
UP = UP[,c("logFC" , "name")]
DOWN = DOWN[,c("logFC" , "name")]

UP$Group = 1
DOWN$Group = 2
all = rbind(UP,DOWN)
UP2 = UP
UP2$logFC = abs(UP2$logFC)
all3 = rbind(UP2,DOWN)


col <- c(rep("red",40), rep("blue",41))
pdf("~/Desktop/newWork/slide/updown1.pdf",width=10,height=10)
ggplot(all , aes( x = logFC  )) + 
  #geom_histogram(binwidth = 0.3 , fill="#DC143C" , colour = "black") +
  geom_histogram(binwidth = 0.3 ,data=subset(all,Group==1), fill="#DC143C" , colour = "black") +
  geom_histogram(binwidth = 0.3 ,data=subset(all,Group==2), fill="#1E90FF" , colour = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.border = element_rect(fill = NA, colour = "black", size=2),
        strip.text.x = element_text(size = 25, colour = "black"),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5)
  ) + xlab("Log Fold Change") +
  ylab("Number of miRNAs") +
  scale_x_continuous(breaks=-5:5) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14)) +
  xlim( c(-5,5))
dev.off()




pdf("~/Desktop/newWork/slide/updown2.pdf",width=10,height=10)
ggplot(all3 , aes( x = logFC  )) + 
  geom_histogram(binwidth = 0.3 , fill="#363636" , colour = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.border = element_rect(fill = NA, colour = "black", size=2),
        strip.text.x = element_text(size = 25, colour = "black"),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_text(size = 20 * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
        axis.text.y =       element_text(size = 20 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x = element_text(size = 20 * 0.8, lineheight = 0, colour = "black", vjust = 0),
        axis.title.y = element_text(size = 20 * 0.8, lineheight = 1, colour = "black", hjust = 0.5)
  ) + xlab("Log Fold Change") +
  ylab("Number of miRNAs") +
  scale_x_continuous(breaks=-5:5) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24)) +
  xlim( c(1,5))
dev.off()




################# historgram for averange 

all_selected
all_selected$name = row.names(all_selected)

pdf("~/Desktop/newWork/slide/average.pdf",width=17,height=10)
ggplot(all_selected , aes( x = name , y = mean )) + 
  geom_point(stat = "identity",size = 8) +
  scale_x_discrete(limits=row.names(all_selected) )+
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.border = element_rect(fill = NA, colour = "black", size=2),
        strip.text.x = element_blank(),
        plot.background = element_rect(fill = '#FFE4E1'),
        strip.background = element_rect(fill = "#FFE4E1"),
        axis.text.x =       element_blank(),
        axis.text.y =       element_text(size = 40 * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
        axis.title.x =  element_text(size = 40 * 0.8, lineheight = 0.9, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 40 * 0.8, lineheight = 1, colour = "black", hjust = 0.5),
        legend.position = "none"
  ) + xlab("Differential Expression Candidates") +
  ylab("Average Expression") 

dev.off()

