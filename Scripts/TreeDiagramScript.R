dbasel_DE_function<-function(Control_AMvsPM, Cvsh_AM, Cvsh_PM){
  require(dplyr)
  controlAMvsPM<-Control_AMvsPM[which(Control_AMvsPM[,5]<0.05),]
  HeatAM<-Cvsh_AM[which(Cvsh_AM[,5]<0.05),]
  HeatPM<-Cvsh_PM[which(Cvsh_PM[,5]<0.05),]
  
  HeatAM<-cbind(HeatAM, rownames(HeatAM))
  HeatPM<-cbind(HeatPM, rownames(HeatPM))
  controlAMvsPM<-cbind(controlAMvsPM, rownames(controlAMvsPM))
  colnames(HeatPM)<-c("logFC", "logCPM", "F", "PValue","FDR", "gene")
  colnames(HeatAM)<-c("logFC", "logCPM", "F", "PValue","FDR", "gene")
  colnames(controlAMvsPM)<-c("logFC", "logCPM", "F", "PValue","FDR", "gene")
  
  HeatBoth<-inner_join(HeatAM, HeatPM, by="gene")
  #print(head(HeatBoth))
  uAM_H<-anti_join(HeatAM, HeatBoth, by="gene")
  uPM_H<-anti_join(HeatPM, HeatBoth, by="gene")
  
  dbasel_higherAM<-subset(controlAMvsPM, logFC < -.5)
  #print(head(dbasel_higherAM))
  Dbasel_higherAM_DE_AM<-semi_join(uAM_H, dbasel_higherAM, by ="gene")
  Dbasel_higherAM_DE_upAM<-subset(Dbasel_higherAM_DE_AM, logFC > 0.5)
  Dbasel_higherAM_DE_downAM<-subset(Dbasel_higherAM_DE_AM, logFC < -0.5)
  
  Dbasel_higherAM_DE_PM<-semi_join(uPM_H, dbasel_higherAM, by ="gene")
  Dbasel_higherAM_DE_upPM<-subset(Dbasel_higherAM_DE_PM, logFC > 0.5)
  Dbasel_higherAM_DE_downPM<-subset(Dbasel_higherAM_DE_PM, logFC < -0.5)
  
  dbasel_higherPM<-subset(controlAMvsPM, logFC > .5)
  #print(head(dbasel_higherPM))
  Dbasel_higherPM_DE_AM<-semi_join(uAM_H, dbasel_higherPM, by ="gene")
  Dbasel_higherPM_DE_upAM<-subset(Dbasel_higherPM_DE_AM, logFC > 0.5)
  Dbasel_higherPM_DE_downAM<-subset(Dbasel_higherPM_DE_AM, logFC < -0.5)
  
  Dbasel_higherPM_DE_PM<-semi_join(uPM_H, dbasel_higherPM, by ="gene")
  Dbasel_higherPM_DE_upPM<-subset(Dbasel_higherPM_DE_PM, logFC > 0.5)
  Dbasel_higherPM_DE_downPM<-subset(Dbasel_higherPM_DE_PM, logFC < -0.5)
  
  nobasel_DE_AM<-anti_join(uAM_H, controlAMvsPM, by="gene")
  nobasel_DE_upAM<-subset(nobasel_DE_AM, logFC > 0.5)
  nobasel_DE_downAM<-subset(nobasel_DE_AM, logFC < -0.5)
  
  nobasel_DE_PM<-anti_join(uPM_H, controlAMvsPM, by="gene")
  nobasel_DE_upPM<-subset(nobasel_DE_PM, logFC > 0.5)
  nobasel_DE_downPM<-subset(nobasel_DE_PM, logFC < -0.5)
  
  ###heat both .x is AM and .y is PM
  
  dbasel_higherAM_DE_BOTH<-left_join(HeatBoth, dbasel_higherAM, by="gene")
  dbasel_higherAM_DE_BOTH<-dbasel_higherAM_DE_BOTH[complete.cases(dbasel_higherAM_DE_BOTH),]
  #print(head(dbasel_higherAM_DE_BOTH))
  dbasel_higherPM_DE_BOTH<-left_join(HeatBoth, dbasel_higherPM, by="gene")
  dbasel_higherPM_DE_BOTH<-dbasel_higherPM_DE_BOTH[complete.cases(dbasel_higherPM_DE_BOTH),]
  #print(head(dbasel_higherPM_DE_BOTH))
  nobasel_DE_BOTH<-anti_join(HeatBoth, controlAMvsPM, by="gene")
  
  Dbasel_higherAM_DE_upBOTH<-subset(dbasel_higherAM_DE_BOTH, logFC.x > 0.5 & logFC.y > 0.5)
  Dbasel_higherPM_DE_upBOTH<-subset(dbasel_higherPM_DE_BOTH, logFC.x > 0.5 & logFC.y > 0.5)
  nobasel_DE_upBOTH<-subset(nobasel_DE_BOTH, logFC.x > 0.5 & logFC.y > 0.5)
  
  Dbasel_higherAM_DE_downBOTH<-subset(dbasel_higherAM_DE_BOTH, logFC.x < -0.5 & logFC.y < -0.5)
  Dbasel_higherPM_DE_downBOTH<-subset(dbasel_higherPM_DE_BOTH, logFC.x < -0.5 & logFC.y < -0.5)
  nobasel_DE_downBOTH<-subset(nobasel_DE_BOTH, logFC.x < -0.5 & logFC.y < -0.5)
  
  Dbasel_higherAM_DE_upAMdownPM<-subset(dbasel_higherAM_DE_BOTH, logFC.x > 0.5 & logFC.y < -0.5)
  Dbasel_higherPM_DE_upAMdownPM<-subset(dbasel_higherPM_DE_BOTH, logFC.x > 0.5 & logFC.y < -0.5)
  nobasel_DE_upAMdownPM<-subset(nobasel_DE_BOTH, logFC.x > 0.5 & logFC.y < -0.5)
  
  Dbasel_higherAM_DE_downAMdupPM<-subset(dbasel_higherAM_DE_BOTH, logFC.x < -0.5 & logFC.y > 0.5)
  Dbasel_higherPM_DE_downAMdupPM<-subset(dbasel_higherPM_DE_BOTH, logFC.x < -0.5 & logFC.y > 0.5)
  nobasel_DE_downAMdupPM<-subset(nobasel_DE_BOTH, logFC.x < -0.5 & logFC.y > 0.5)
  
  group_list<-list(hAM_DE_upAM=Dbasel_higherAM_DE_upAM, hAM_DE_dwAM=Dbasel_higherAM_DE_downAM, hAM_DE_upPM=Dbasel_higherAM_DE_upPM, hAM_DE_dwPM=Dbasel_higherAM_DE_downPM, hAM_DE_upBOTH=Dbasel_higherAM_DE_upBOTH, hAM_DE_dwBOTH=Dbasel_higherAM_DE_downBOTH, hAM_DE_upAMdwPM=Dbasel_higherAM_DE_upAMdownPM, hAM_DE_dwAMupPM=Dbasel_higherAM_DE_downAMdupPM,
                   hPM_DE_upAM=Dbasel_higherPM_DE_upAM, hPM_DE_dwAM=Dbasel_higherPM_DE_downAM, hPM_DE_upPM=Dbasel_higherPM_DE_upPM, hPM_DE_dwPM=Dbasel_higherPM_DE_downPM, hPM_DE_upBOTH=Dbasel_higherPM_DE_upBOTH, hPM_DE_dwBOTH=Dbasel_higherPM_DE_downBOTH, hPM_DE_upAMdwPM=Dbasel_higherPM_DE_upAMdownPM, hPM_DE_dwAMupPM=Dbasel_higherPM_DE_downAMdupPM,
                   eq_DE_upAM=nobasel_DE_upAM, eq_DE_dwAM=nobasel_DE_downAM, eq_DE_upPM=nobasel_DE_upPM, eq_DE_dwPM=nobasel_DE_downPM, eq_DE_upBOTH=nobasel_DE_upBOTH, eq_DE_dwBOTH=nobasel_DE_downBOTH, eq_DE_upAMdwPM=nobasel_DE_upAMdownPM, eq_DE_dwAMupPM=nobasel_DE_downAMdupPM)
  for(i in 1:length(group_list)){
    n<-names(group_list)
    print(n[i])
    print(length(group_list[[i]][,1]))
  }
  return(group_list)
}