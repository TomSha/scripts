prefix_list<-read.table("prefix_list_rep1")$V1

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
	if(nEP==34){
		epochs_list[[i]]<-ep[,2][2:nEP*2]
	}else{
		epochs_list[[i]]<-ep[,2][1:nEP*2]
	}
}

index_sizes_list<-lapply(epochs_list, grep, pattern="SF15")

output_sizes_list<-mapply(function(output,index) output[,c(1,index+1)], outputRCT_list, index_sizes_list, SIMPLIFY=F)

output_sizes_list<-lapply(output_sizes_list,function(x) t(apply(x,1,function(x) c(x[1],x[2:length(x)]/mean(x[2:length(x)])))))

epochs_sizes_list<-lapply(epochs_list,grep,pattern="SF15",value=T)

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i]
	folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="");
	folder_sizes=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox_DS/",sep="");

	if (!dir.exists(folder_sizes)) dir.create(folder_sizes)
	
	file.copy(paste(folder,prefix,".bin",sep=""),paste(folder_sizes,prefix,".bin",sep=""))
	print(paste("copied binary",prefix))
	file.copy(paste(folder,prefix,"_cell_centers.dat",sep=""),paste(folder_sizes,prefix,"_cell_centers.dat",sep=""))

	file.copy(paste(folder,prefix,"_cleaned.nii",sep=""),paste(folder_sizes,prefix,"_cleaned.nii",sep=""))

	file.copy(paste(folder,prefix,"_noise_thresh.dat",sep=""),paste(folder_sizes,prefix,"_noise_thresh.dat",sep=""))

	file.copy(paste(folder,prefix,"_mip.dat",sep=""),paste(folder_sizes,prefix,"_mip.dat",sep=""))

	file.copy(paste(folder,prefix,"_corvec.dat",sep=""),paste(folder_sizes,prefix,"_corvec.dat",sep=""))

	file.copy(paste(folder,prefix,"_output_rawCT.dat",sep=""),paste(folder_sizes,prefix,"_output_rawCT.dat",sep=""))

	write.table(x=output_sizes_list[[i]],file=paste(folder_sizes,prefix,"_outputCTM.dat",sep=""),row.names=F,col.names=F)
	write.table(x=epochs_sizes_list[[i]],file=paste(folder_sizes,prefix,"_epochsC.dat",sep=""),row.names=F,col.names=F)
	

	print(paste("finished",prefix))
}




