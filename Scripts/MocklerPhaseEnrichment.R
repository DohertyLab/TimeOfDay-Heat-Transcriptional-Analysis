MocklerEnrichReg<-function(NAME,CUT)
{
  MocklerFull$V2<-gsub(",","",MocklerFull$V2)
  MocklerSubset<-MocklerFull[which(MocklerFull$V4==NAME),]
  MocklerSubset<-MocklerSubset[which(MocklerSubset$V2!=""),]
  MocklerCut<-MocklerSubset[which(MocklerSubset$V6>CUT),]
  Overlap<-intersect(MocklerCut$V2,HSG_Regulators)
  x1<-length(Overlap)
  m1<-length(MocklerCut$V1)
  n1<-length(MocklerSubset$V1)-length(MocklerCut$V1)
  k1<-length(HSG_Regulators)
  OMG<-phyper(x1,m1,n1,k1, lower.tail = FALSE)
  return(OMG)
}

MocklerEnrichResp<-function(NAME,CUT)
{
  MocklerFull$V2<-gsub(",","",MocklerFull$V2)
  MocklerSubset<-MocklerFull[which(MocklerFull$V4==NAME),]
  MocklerSubset<-MocklerSubset[which(MocklerSubset$V2!=""),]
  MocklerCut<-MocklerSubset[which(MocklerSubset$V6>CUT),]
  Overlap<-intersect(MocklerCut$V2,HSG_Responders)
  x1<-length(Overlap)
  m1<-length(MocklerCut$V1)
  n1<-length(MocklerSubset$V1)-length(MocklerCut$V1)
  k1<-length(HSG_Responders)
  OMG<-phyper(x1,m1,n1,k1, lower.tail = FALSE)
  return(OMG)
}