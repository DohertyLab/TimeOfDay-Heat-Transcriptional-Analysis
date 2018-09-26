Grepper<-function(GO_Number,Data)
{
  Result<-Data[grep(GO_Number,names(Data))]
  print(Data[1])
  print(Result)
}