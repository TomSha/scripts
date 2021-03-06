args <- commandArgs(trailingOnly = TRUE)

prefix_list_file = args[1]
folder_pulled = args[2]



require(plotrix)
require(mmand)

#mainresultfolder="/home/meyer-lab/SeG/results_toolbox/"
prefix_list<-as.character(read.table(prefix_list_file)$V1);
centers_list<-vector("list",length(prefix_list));
centers_sd_list<-vector("list",length(prefix_list));
mainfolder="/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/"

resc <- function(x,zero=0.2){
	return(zero+(0.999-zero)*(x-min(x))/(max(x)-min(x)))
}

for(pid in 1:length(prefix_list)){
	prefix=prefix_list[pid]
	folder=paste(mainfolder,prefix,"/results_toolbox/",sep="")
#	folder=paste(mainfolder,prefix,"/results_toolbox/",sep="")
	epfile=paste(folder,prefix,"_rslice10_epochsC.dat",sep="");
	nifti=paste(mainfolder,"/",prefix,"/nii/rslice10/rim10.nii",sep="");
	
	#how many nearest neighbours need to be in same cluster for cell not to be removed as a halo
	KNN=2
	EP<-read.table(epfile);
	
	nEP=nrow(EP)
	EPlist<-vector(mode="character", length=nEP);
	EPlist<-EP$V1

	data0=read.table(paste(folder,prefix,"_combined_outputC.dat",sep=""));
	dencl<-read.table(paste(folder,prefix,"_summary.dat",sep=""));
	KNNR0<-read.table(paste(folder,prefix,"_KNNrecords.dat",sep=""))[,1:(KNN+1)];
	thsum<-as.numeric(read.table(paste(folder,prefix,"_thresholds.dat",sep="")));
	cenfile<-read.table(paste(folder,prefix,"_cell_centers.dat",sep=""));
	
	cenfile<-cenfile[data0$V2 %in% dencl$V1,]
	data=data0[data0$V2 %in% dencl$V1,];
	#V1 cells V2+ KNNs
	KNNR=KNNR0[KNNR0$V1 %in% dencl$V1,];
	
	dencl <- dencl[order(dencl$V1),];
	cl<-dencl$V6+1;
	
	#creates list of cluster membership of a cell's KNNs
	halomatrix<-matrix(ncol=KNN,nrow=nrow(dencl));
	for(i in 1:nrow(dencl)){
		halomatrix[i,]<-cl[as.numeric(KNNR[i,2:(KNN+1)]+1)]
	}
	
	hm<-cbind(halomatrix,cl)
	
	#if KNNs not in same cluster that cell == 0
	hl<-(apply(hm,1,function (x) max(x)-min(x)))==0
	listminprob<-vector("numeric",max(cl))
	
	#Take cells in halo (!hl) and find most dense cell, this will be min density threshold for which cells must have to be kept in cluster
	for(i in 1:max(cl)) {
		if(nrow(dencl[cl==i & !hl,])>0) {
			listminprob[i]=max(dencl$V3[cl==i & (!hl)]);
		} else {
			listminprob[i]=0;
		}
	}
	threshold<-(hl & dencl$V3>listminprob[cl])
	write.table(file=paste(folder,prefix,"_inc.dat",sep=""),threshold)

	
	sizes<-rep(0,max(cl));
	sizes[as.numeric(names(table(cl[threshold])))]=as.numeric(table(cl[threshold]));sizes_cl<-sizes;
	
	#centers are the average response of each epoch for each cluster
	centers<-matrix(nrow=max(cl),ncol=nEP);
	centers_sd<-matrix(nrow=max(cl),ncol=nEP);
	
	#create data frame of output.dat for cells in a given cluster above threshold, create mean and sd of cells in cluster
	for(i in 1:max(cl)){
		tmpframe=data[cl==i & threshold,1:nEP+2];
		tmpdencl=dencl[cl==i & threshold,];
		if(!is.null(nrow(tmpframe))){
			if(nrow(tmpframe)==0){
				centers[i,] <- rep(NA,nEP);
				centers_sd[i,] <- rep(NA,nEP);
			}
			if(nrow(tmpframe)>1){
				centers[i,] <- as.numeric(tmpframe[which.max(tmpdencl$V3),])
				centers_sd[i,] <- apply(tmpframe,2,sd);
			}
		} else {
			centers[i,] <- rep(NA,nEP);
			centers_sd[i,] <- rep(NA,nEP);
		}
	}
	
	
	colnames(centers)<-EPlist
	colnames(centers_sd)<-EPlist
	write.table(x=data.frame(PREFIX=prefix,SIZE=sizes,centers),file=paste(folder,prefix,"_centers.dat",sep=""),row.names=F)
	write.table(x=data.frame(PREFIX=prefix,SIZE=sizes,centers_sd),file=paste(folder,"/",prefix,"_centers_sd.dat",sep=""),row.names=F)
	centers_list[[pid]]=data.frame(PREFIX=prefix,SIZE=sizes,CLID=1:max(cl),centers)
	centers_sd_list[[pid]]=data.frame(PREFIX=prefix,SIZE=sizes,CLID=1:max(cl),centers_sd)
	selec=complete.cases(centers_list[[pid]]);
	centers_list[[pid]]=centers_list[[pid]][selec,];
	centers_sd_list[[pid]]=centers_sd_list[[pid]][selec,];

}

all_centers<-do.call(rbind,centers_list)
all_centers_sd<-do.call(rbind,centers_sd_list)

write.table(x=data.frame(ID=1:nrow(all_centers)-1,all_centers),file=paste(folder_pulled,"/all_centers.dat",sep=""),row.names=F,col.names=F)
write.table(x=data.frame(ID=1:nrow(all_centers_sd)-1,all_centers_sd),file=paste(folder_pulled,"/all_centers_sd.dat",sep=""),row.names=F,col.names=F)
write.table(x=data.frame(colnames(all_centers)[3+1:nEP]),file=paste(folder_pulled,"/pulled_epochsC.dat",sep=""),col.names=F,row.names=F)
#write.table(x=data.frame(colnames(all_centers)[3+1:nEP]),file=paste(folder_pulled,"/pulled_epochs.dat",sep=""),col.names=F,row.names=F)
