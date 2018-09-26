GOMasterPlotter<-function(GOTerm)
{
  NumeroMaximo<-max(sapply(BIGLISTCLEAN,function (y) -log10(y[GOTerm])))
  barplot(sapply(BIGLISTCLEAN,function (y) -log10(y[GOTerm])),las=2,names.arg = names(BIGLISTCLEAN),ylim=c(0,NumeroMaximo+1),main = AllGONames[GOTerm],cex.names=0.8)
  abline(h=1.3,col="red",lty=2)
}

