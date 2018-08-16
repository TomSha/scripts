args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
resultsfolder <- args[2]

require(plotrix)
require(mmand)

output=read.table(paste(resultsfolder,"/",prefix,"_output.dat",sep=""));
epochs=read.table(paste(args[3],prefix,"_time.log",sep=""))$V2;
timelog=read.table(paste(args[3],prefix,"_time.log",sep=""))$V1+90;

#---------remove epochs from raw output---------------
outputR=read.table(paste(resultsfolder,"/",prefix,"_output_raw.dat",sep=""));
#------------------------------------------------------

nEP=length(epochs)/2;
epochs=epochs[1:nEP *2];

outputfolderC <- paste(resultsfolder,"/",prefix,"_outputC.dat",sep="");
epochfolderC <- paste(resultsfolder,"/",prefix,"_epochsC.dat",sep="");
outputfolderRC <- paste(resultsfolder,"/",prefix,"_output_rawC.dat",sep="");
EPtimefolderC <- paste(resultsfolder,"/",prefix,"_EPtimeC.dat",sep="");
temp <- paste(resultsfolder,"/",prefix,"_temp.dat",sep="");


#---------------------------remove epoch responses------------------------------
stim<-"WIGGLY|CONCENTRIC|LOOM"
epochsI<- grep((stim), epochs, value = F,invert=T);
outputC<-output[,c(1,epochsI+1)];
#cmean <- colMeans(outputC[,2:ncol(outputC)]);
#outputC[,2:ncol(outputC)] <- outputC[,2:ncol(outputC)]/cmean[col(outputC[,2:ncol(outputC)])];
outputC[,2:ncol(outputC)] <- t(apply(outputC[,2:ncol(outputC)],1,function(x) (x-mean(x))/sd(x))) 

write.table(outputC, outputfolderC, row.names = F, col.names = F);


outputRC<-outputR[,c(1,epochsI+1)];
write.table(outputRC, outputfolderRC, row.names = F, col.names = F);

#--------------Generate vector of names of epochs----------------------
epochsC<- grep((stim), epochs, value = T,invert=T);
write.table(epochsC, epochfolderC, row.names = F, col.names = F);

#---calculate epoch frame number---
timelogM=matrix(floor(timelog*20),ncol=2,byrow=T);
EPtime=timelogM[,2]-timelogM[,1];
EPtimeC=EPtime[epochsI];
write.table(EPtimeC,EPtimefolderC,row.names=F,col.names=F);

