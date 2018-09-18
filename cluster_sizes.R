prefix_list<-read.table("prefix_list_bars_dots")$V1

folder_list<-vector("list",length(prefix_list))
outputRCT_list<-vector("list",length(prefix_list))
epochs_list<-vector("list",length(prefix_list))

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i]
	folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="");
	epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="");
	
	folder_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",prefix,sep="");
	outputRCT_list[[i]]<-read.table(paste(folder,prefix,"_output_rawCT.dat",sep=""));
	ep<-read.table(epfile);
	nEP<-length(ep$V2)/2
	epochs_list[[i]]<-ep[,2][2:nEP*2]
}

index_dot_90_list<-lapply(epochs_list, grep, pattern="DOT_90")
index_dot_270_list<-lapply(epochs_list, grep, pattern="DOT_270")

output_dot_90_list<-mapply(function(output,index) output[,c(1,index+1)], outputRCT_list, index_dot_90_list, SIMPLIFY=F)
output_dot_270_list<-mapply(function(output,index) output[,c(1,index+1)], outputRCT_list, index_dot_270_list, SIMPLIFY=F)

output_dot_90_list<-lapply(output_dot_90_list,function(x) t(apply(x,1,function(x) c(x[1],x[2:length(x)]/mean(x[2:length(x)])))))
output_dot_270_list<-lapply(output_dot_270_list,function(x) t(apply(x,1,function(x) c(x[1],x[2:length(x)]/mean(x[2:length(x)])))))

epochs_dot_90_list<-lapply(epochs_list,grep,pattern="DOT_90",value=T)
epochs_dot_270_list<-lapply(epochs_list,grep,pattern="DOT_270",value=T)

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i]
	folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="");

	write.table(x=output_dot_90_list[[i]],file=paste(folder,prefix,"_output_dot_90CTM.dat",sep=""),row.names=F,col.names=F)
	write.table(x=output_dot_270_list[[i]],file=paste(folder,prefix,"_output_dot_270CTM.dat",sep=""),row.names=F,col.names=F)
	write.table(x=epochs_dot_90_list[[i]],file=paste(folder,prefix,"_epochs_dot_90.dat",sep=""),row.names=F,col.names=F)
	write.table(x=epochs_dot_270_list[[i]],file=paste(folder,prefix,"_epochs_dot_270.dat",sep=""),row.names=F,col.names=F)
}

