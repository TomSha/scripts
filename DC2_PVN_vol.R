#args <- commandArgs(trailingOnly = TRUE)

#prefix<- args[1]


#CHANGED 1  060818 removed min density for halo removal

require(plotrix)
require(mmand)

prefix="171012_PVN_barrage_F1_1_v10"
#folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox_sizes/",sep="")
folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox/",sep="")
epfile=paste(folder,prefix,"_rslice10_epochsC.dat",sep="");
nifti=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/nii/rslice10/rim10.nii",sep="");

#how many nearest neighbours need to be in same cluster for cell not to be removed as a halo
KNN=1
EP<-read.table(epfile);
nEP<-nrow(EP)
EPlist<-vector(mode="character", length=nEP);
EPlist<-EP$V1

#data0=read.table(paste(folder,prefix,"_combined_output_dot_270CTM.dat",sep=""));
#noisethresh<-read.table(paste(folder,prefix,"_noise_thresh.dat",sep=""))$V1;
#dencl<-read.table(paste(folder,prefix,"_summary_dot_270.dat",sep=""));
#KNNR0<-read.table(paste(folder,prefix,"_KNNrecords_dot_270.dat",sep=""))[,1:(KNN+1)];
#thsum<-as.numeric(read.table(paste(folder,prefix,"_thresholds_dot_270.dat",sep="")));

data0=read.table(paste(folder,prefix,"_combined_outputC.dat",sep=""));
noisethresh<-read.table(paste(folder,prefix,"_noise_thresh.dat",sep=""))$V1;
dencl<-read.table(paste(folder,prefix,"_summary.dat",sep=""));
KNNR0<-read.table(paste(folder,prefix,"_KNNrecords.dat",sep=""))[,1:(KNN+1)];
thsum<-as.numeric(read.table(paste(folder,prefix,"_thresholds.dat",sep="")));
cenfile0<-read.table(paste(folder,prefix,"_cell_centers.dat",sep=""))[noisethresh,];

cenfile<-cenfile0[data0$V2 %in% dencl$V1,]
data=data0[data0$V2 %in% dencl$V1,];
#V1 cells V2+ KNNs
KNNR=KNNR0[KNNR0$V1 %in% dencl$V1,];

#cl=cluster membership
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
#threshold - cells with KNNs in same cluster and have density higher than threshold
#threshold<-(hl & dencl$V3>listminprob[cl])#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHANGED 1
threshold<-(hl)
#ord<-order(EPlist)
ord<-order(as.numeric(gsub("\\D","",EPlist)))

sizes<-rep(0,max(cl));
sizes[as.numeric(names(table(cl[threshold])))]=as.numeric(table(cl[threshold]));sizes_cl<-sizes;
ordcl<-order(sizes,decreasing=T);

#centers are the average response of each epoch for each cluster
centers<-matrix(nrow=max(cl),ncol=nEP);
centers_sd<-matrix(nrow=max(cl),ncol=nEP);
colors=rainbow(max(cl));

#create data frame of output.dat for cells in a given cluster above threshold, create mean and sd of cells in cluster
for(i in 1:max(cl)){
	tmpframe=data[cl==i & threshold,1:nEP+2];
	tmpdencl=dencl[cl==i & threshold,];
	if(!is.null(nrow(tmpframe))){
		if(nrow(tmpframe)==0){
			centers[i,] <- rep(0,nEP);
			centers_sd[i,] <- rep(0,nEP);
		}
		if(nrow(tmpframe)>1){
			centers[i,] <- as.numeric(tmpframe[which.max(tmpdencl$V3),])
			centers_sd[i,] <- apply(tmpframe,2,sd);
		}
	} else {
		centers[i,] <- tmpframe;
		centers_sd[i,] <- rep(0,nEP);
	}
}


cell_id<-function(rowNum){
	rowNum<-as.numeric(rowNum)
	slice<-data[rowNum,1]
	if (!slice == min(data[,1])){
		subtract<-max(which(data[,1]==slice-1))
		cellID<-data[rowNum,2]-data[subtract,2]
	}else{
		cellID<-data[rowNum,2]+1
	}
	returnList<-list("slice"=slice,"cellID"=cellID)
	return(returnList)
}



plot_dp <- function(){
	plot(log(dencl$V4)[threshold],log(dencl$V3)[threshold],col=colors[ordcl[cl[threshold]]],pch=19,cex=0.6,xlab="log(distance)",ylab="log(density)");
	lines(rep(thsum[1],100),seq(-10,10,length=100),lty=2,col="green")
	lines(seq(-10,10,length=100),-nEP*seq(-10,10,length=100)+thsum[2],lty=2,col="green")
	lines(seq(-100,100,length=100),rep(thsum[3],100),lty=2,col="green")
}

plot_dp_free <- function(c){
	logdist<-seq(min(log(dencl$V4)),max(log(dencl$V4)),length.out=100);
	logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
	colors<-rep("black",nrow(dencl));
	colors[cl==c]="red";
	plot(log(dencl$V4),log(dencl$V3),col=colors,pch=18);
	lines(rep(thsum[1],100),logprob,lty=2,col="green")
	lines(seq(-100,100,length=100),rep(thsum[3],100),lty=2,col="green")
	lines(logdist,-nEP*logdist+thsum[2],lty=2,col="green")
}

plot_traj<-function(X,min=0,max=1,NS=1,AS=1){
	par(las=2,mar=c(3,10,2,3))
	ord2<-order(centers[X,]);
	pos<-barplot(centers[X,ord2],names.arg=EPlist[ord2],horiz=T,xlim=c(0,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS,main=X);
	arrows(y0=pos,y1=pos,x0=centers[X,ord2]-centers_sd[X,ord2],x1=centers[X,ord2]+centers_sd[X,ord2],code=3,length=0);

}

plot_box<-function(X){
	par(mar=c(3,10,2,2))
#	ord2<-order(centers[X,]);
#	boxplot(data[cl==X & threshold,2+ord2],names=EPlist[ord2],las=2,horizontal=T,main=paste("cluster",X,":" ,sizes[X],"cells"))
	boxplot(data[cl==X & threshold,2+ord],names=EPlist[ord],las=2,horizontal=T,main=paste("cluster",X,":" ,sizes[X],"cells"))
#	abline(v=0, lty=2,col="blue")
}
write.table(data.frame(cenfile[threshold,],cl[threshold],dencl$V3[threshold]),file=paste(folder,prefix,"_plot3d.dat",sep=""),col.names=F,row.names=F);

plot_cell<-function(id,xmax=100000){
	ep<-read.table(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep=""))
	ep$V1<-(ep$V1+90)*4
	nEP2<-nrow(ep)/2
	param<-cell_id(id);
	print(param);
	slice<-paste("_rslice",param$slice,sep="");
	par(mar=c(8,3,2,1));
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",param$cellID," ",paste(folder,prefix,sep=""),slice,sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix,id),xaxt='n',ann=F);
	axis(1,at=ep$V1[1:nEP2*2-1],label=ep$V2[1:nEP2*2-1],las=2);
	rect(xleft=ep$V1[1:nEP2 *2 -1],xright=ep$V1[1:nEP2 *2 ],ybottom=rep(-1000,nEP2*2),ytop=rep(xmax,nEP2*2),col=rgb(0,1,0,.05));

}

plot_cell_in_cluster<-function(id,cluster){
	lmat<-matrix(c(1,2,3,3),ncol=2,byrow=T);
	layout(lmat);
	plot_box(cluster);
	par(mar=c(5,2,2,2))
	barcent<-barplot(as.numeric(data0[threshold,][cl[threshold]==cluster,][id,3:ncol(data0)]))
	axis(1,at=barcent,labels=EPlist,las=2,cex.axis=0.8)
	cell_row<-which(data0[,2]==data0[cl==cluster,][id,2])
	plot_cell(cell_row)

}

plot_fig<-function(mat=c(1,2,3,4,5,6,7),cluster,cols=1){
	lmat<-matrix(mat,ncol=cols,byrow=T)
	layout(lmat);
	par(mar=c(5,2,2,2));
	plot_box(cluster);
	for (id in 1:sizes[cluster]){
		 cell_row<-which(data0[,2]==data0[cl==cluster,][id,2])
		 plot_cell(cell_row)
		 }
}


