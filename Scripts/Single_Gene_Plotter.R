GGG_PlotDaddy<-function(Gene,Y_Axis_Scaler,Gene_Name)
{
  GeneRow<-CPM_Table[Gene,]
  Means<-c(mean(GeneRow[1:4]),mean(GeneRow[9:12]),mean(GeneRow[5:8]),mean(GeneRow[13:16]))
  SDs<-c(sd(GeneRow[1:4])/sqrt(length(GeneRow[1:4])),sd(GeneRow[9:12])/sqrt(length(GeneRow[9:12])),sd(GeneRow[5:8])/sqrt(length(GeneRow[5:8])),sd(GeneRow[13:16])/sqrt(length(GeneRow[13:16])))
  Treatment_Labels<-c("C_AM","H_AM","C_PM","H_PM")
  PlotReady<-cbind(Means,SDs,Treatment_Labels)
  LetsPlot<-as.data.frame(PlotReady)
  LetsPlot$Means <- as.numeric(as.character(LetsPlot$Means))
  LetsPlot$SDs <- as.numeric(as.character(LetsPlot$SDs))
  LetsPlot$Treatment_Labels <- factor(LetsPlot$Treatment_Labels, levels = c("C_AM","H_AM","C_PM","H_PM"))
  
  ggplot(data=LetsPlot, aes(x=Treatment_Labels, y=Means, fill=Treatment_Labels)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_bw()+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),axis.text.y = element_text(face="bold", size=22),plot.title = element_text(hjust = 0.5,size = 60,face = "bold"),panel.grid.major.y=element_line(color="darkgrey"))+
    scale_fill_manual(values=c('firebrick1','darkred','dodgerblue','dodgerblue4'),labels=c("Dawn Control","Dawn Heat","Dusk Control","Dusk Heat"),name="Treatment")+
    geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs),
                  width=.2,
                  position=position_dodge(0.9))+
    scale_y_continuous(expand=c(0,0))+
    expand_limits(y=Y_Axis_Scaler)+
    guides(fill=FALSE)+
    ggtitle(Gene_Name)
}