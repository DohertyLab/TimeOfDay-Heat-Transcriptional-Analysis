#inputs list of pval vectors and pulls out GO term of values with valid go term names
CleanPVALVectors<-function(pvals,RegularExp)
{
  GONames<-str_extract(string = names(pvals),pattern = RegularExp)
  names(pvals)<-GONames
  Keep<-GONames[!is.na(GONames)]
  pvals<-pvals[Keep]
  pvals<-as.numeric(pvals)
  names(pvals)=Keep
  return(pvals)
}