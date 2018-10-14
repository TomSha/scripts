args <- commandArgs(trailingOnly=TRUE)

folder = args[1]
gprefix = args[2]
init_slice=as.numeric(args[3])
final_slice=as.numeric(args[4])
data_list <- vector("list",final_slice-init_slice+1)
data_list_raw <- vector("list",final_slice-init_slice+1)
data_list_max <- vector("list",final_slice-init_slice+1)
cell_centers<-vector("list",final_slice-init_slice+1)
noise_thresh<-vector("list",final_slice-init_slice+1)
for(i in init_slice:final_slice){
	cat(i,"\n")
	data_list[[i]] <- read.table(paste(folder,"/",gprefix,"_rslice",i,"_outputCT.dat",sep=""))
	data_list_raw[[i]] <- read.table(paste(folder,"/",gprefix,"_rslice",i,"_output_rawCT.dat",sep=""))
	data_list_max[[i]] <- read.table(paste(folder,"/",gprefix,"_rslice",i,"_output_maxCT.dat",sep=""))
	cell_centers[[i]] <- read.table(paste(folder,"/",gprefix,"_rslice",i,"_cell_centers.dat",sep=""))
	noise_thresh[[i]]<- read.table(paste(folder,"/",gprefix,"_rslice",i,"_noise_thresh.dat",sep=""))
	data_list[[i]] <- cbind(rep(i,nrow(data_list[[i]])),data_list[[i]])
	data_list_raw[[i]] <- cbind(rep(i,nrow(data_list_raw[[i]])),data_list_raw[[i]])
	data_list_max[[i]] <- cbind(rep(i,nrow(data_list_max[[i]])),data_list_max[[i]])
	cell_centers[[i]] <- cbind(cell_centers[[i]],rep(i,nrow(cell_centers[[i]])))
	noise_thresh[[i]] <- cbind(noise_thresh[[i]],rep(i,nrow(noise_thresh[[i]])))
	if(i>init_slice){
		data_list[[i]][,2]=data_list[[i]][,2]+data_list[[i-1]][nrow(data_list[[i-1]]),2]+1
		data_list_raw[[i]][,2]=data_list_raw[[i]][,2]+data_list_raw[[i-1]][nrow(data_list_raw[[i-1]]),2]+1
		data_list_max[[i]][,2]=data_list_max[[i]][,2]+data_list_max[[i-1]][nrow(data_list_max[[i-1]]),2]+1
		}
}
combined_data<-do.call(rbind,data_list)
combined_data_raw<-do.call(rbind,data_list_raw)
combined_data_max<-do.call(rbind,data_list_max)
combined_cell_centers<-do.call(rbind,cell_centers)
combined_noise_thresh<-do.call(rbind,noise_thresh)

write.table(combined_data,file=paste(folder,"/",gprefix,"_combined_outputCT.dat",sep=""),row.names=F,col.names=F)
write.table(combined_data_raw,file=paste(folder,"/",gprefix,"_combined_output_rawCT.dat",sep=""),row.names=F,col.names=F)
write.table(combined_data_max,file=paste(folder,"/",gprefix,"_combined_output_maxCT.dat",sep=""),row.names=F,col.names=F)
write.table(combined_cell_centers,file=paste(folder,"/",gprefix,"_cell_centers.dat",sep=""),row.names=F,col.names=F)
write.table(combined_noise_thresh,file=paste(folder,"/",gprefix,"_noise_thresh.dat",sep=""),row.names=F,col.names=F)

print(paste("num epochs in slice",ncol(data_list[[i]])))
print(paste("num epochs in combined",ncol(combined_data)))
