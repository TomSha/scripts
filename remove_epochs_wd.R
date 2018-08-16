args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
resultsfolder <- args[2]

require(plotrix)
require(mmand)

output=read.table(paste(resultsfolder,"/",prefix,"_output.dat",sep=""));
epochs=read.table(args[3])$V2;

#---------remove epochs from raw output---------------
outputR=read.table(paste(resultsfolder,"/",prefix,"_output_raw.dat",sep=""));
#------------------------------------------------------

nEP=length(epochs)/2;
epochs=epochs[1:nEP *2];

outputfolderC <- paste(resultsfolder,"/",prefix,"_outputC.dat",sep="");
epochfolderC <- paste(resultsfolder,"/",prefix,"_epochsC.dat",sep="");
outputfolderRC <- paste(resultsfolder,"/",prefix,"_output_rawC.dat",sep="");
temp <- paste(resultsfolder,"/",prefix,"_temp.dat",sep="");

#-------Generate vectors containing the max response of each cell for each WD quadrant---------
#maxVectorDR <-   rep_len(0, nrow(output));
#maxVectorDL <-  rep_len(0, nrow(output));
#maxVectorUR <-  rep_len(0, nrow(output));
#maxVectorUL <- rep_len(0, nrow(output));
#
#
#epochs_DR <-grep(("_DR"), epochs, value = F, invert=F);
#response_DR <- output[,c(1,epochs_DR+1)]
#for (i in 1:nrow(output))
#{
#	  a <- max(response_DR[i,2:4])
#	    maxVectorDR[i] <- a
#}
#
#epochs_DL <-grep(("_DL"), epochs, value = F, invert=F);
#response_DL <- output[,c(1,epochs_DL+1)]
#for (i in 1:nrow(output))
#{
#	  a <- max(response_DL[i,2:4])
#	    maxVectorDL[i] <- a
#}
#
#epochs_UR <- grep(("_UR"), epochs, value = F, invert=F);
#response_UR <- output[,c(1,epochs_UR+1)]
#for (i in 1:nrow(output))
#{
#	  a <- max(response_UR[i,2:4])
#	    maxVectorUR[i] <- a
#}
#
#epochs_UL <- grep(("_UL"), epochs, value = F, invert=F);
#response_UL <- output[,c(1,epochs_UL+1)]
#for (i in 1:nrow(output))
#{
#	      a <- max(response_UL[i,2:4])
#		          maxVectorUL[i] <- a
#}
#


#---------------------------remove epoch responses------------------------------

stim<-"WIGGLY30_"
epochsC<- grep((stim), epochs, value = F,invert=F);
outputC<-output[,c(1,epochsC+1)];
#--------------------
nEP<-length(epochsC);

#outputC3<- cbind(maxVectorDR,maxVectorDL,maxVectorUR,maxVectorUL)
outputC[,1:nEP+1] <- t(apply(outputC[,1:nEP+1],1,function(x) (x-mean(x))/sd(x))) 

write.table(outputC, outputfolderC, row.names = F, col.names = F);


outputRC<-outputR[,c(1,epochsC+1)];
write.table(outputRC, outputfolderRC, row.names = F, col.names = F);

#--------------Generate vector of names of epochs----------------------
epochsC<- grep((stim), epochs, value = T,invert=F);
#epochsC<- c("WIGGLY_DR","WIGGLY_DL","WIGGLY_UR","WIGGLY_UL");
write.table(epochsC, epochfolderC, row.names = F, col.names = F);
