prefix<-"180727_PVN_barrage4_v10_F1_1"

folder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox/",sep="")
epochs<-read.table(paste(folder,prefix,"_rslice10_epochsC.dat",sep=""))$V1

outputRCT<-read.table(paste(folder,prefix,"_combined_output_rawCT.dat",sep=""));

index_dot_90<-grep(pattern="DOT_90",x=epochs)
index_dot_270<-grep(pattern="DOT_270",x=epochs)

output_dot_90<-outputRCT[,c(1,2,index_dot_90+2)]
output_dot_270<-outputRCT[,c(1,2,index_dot_270+2)]

output_dot_90[,3:ncol(output_dot_90_list)]<-output_dot_90_list[,3:ncol(output_dot_90_list)]/rowMeans(output_dot_90_list[,3:ncol(output_dot_90_list)])
output_dot_270[,3:ncol(output_dot_270_list)]<-output_dot_270_list[,3:ncol(output_dot_270_list)]/rowMeans(output_dot_270_list[,3:ncol(output_dot_270_list)])

epochs_dot_90<-epochs[index_dot_90]
epochs_dot_270<-epochs[index_dot_270]

write.table(x=output_dot_90,file=paste(folder,prefix,"_combined_output_dot_90CTM.dat",sep=""),row.names=F,col.names=F)
write.table(x=output_dot_270,file=paste(folder,prefix,"_combined_output_dot_270CTM.dat",sep=""),row.names=F,col.names=F)
write.table(x=epochs_dot_90,file=paste(folder,prefix,"_epochs_dot_90.dat",sep=""),row.names=F,col.names=F)
write.table(x=epochs_dot_270,file=paste(folder,prefix,"_epochs_dot_270.dat",sep=""),row.names=F,col.names=F)

