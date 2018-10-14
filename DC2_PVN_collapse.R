#args <- commandArgs(trailingOnly = TRUE)

#prefix<- args[1]

require(plotrix)
require(mmand)

prefix="180608_PVN_barrage4_F1_1"
folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox_sizes/",sep="")#CHANGES 17/08/2018
#epfilenoC="/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/170623_PVN_barrage_7dpf_F2_1_time.log";
#epfile=paste(folder,prefix,"_epochsC.dat",sep="");
epfile=paste(folder,prefix,"_epochs_dot_270.dat",sep="");
nifti=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");

#how many nearest neighbours need to be in same cluster for cell not to be removed as a halo
KNN=1
EP<-read.table(epfile);
#EPindex<-grep("CONCENTRIC",EP$V2,invert=T);
#EP<-EP[EPindex,]
#nEP<-nrow(EP)/2
nEP<-nrow(EP)
EPlist<-vector(mode="character", length=nEP);
#EPlist<-as.character(EP[seq(1,nrow(EP),2),2])
EPlist<-EP$V1

data0=read.table(paste(folder,prefix,"_outputCTM.dat",sep=""));
#data0=read.table(paste(folder,prefix,"_output.dat",sep=""));
noisethresh<-read.table(paste(folder,prefix,"_noise_thresh.dat",sep=""))$V1;
dencl<-read.table(paste(folder,prefix,"_summary.dat",sep=""));
KNNR0<-read.table(paste(folder,prefix,"_KNNrecords.dat",sep=""))[,1:(KNN+1)];
thsum<-as.numeric(read.table(paste(folder,prefix,"_thresholds.dat",sep="")));
cenfile<-read.table(paste(folder,prefix,"_cell_centers.dat",sep=""))[noisethresh,];
#nii_cleaned<-paste(base_folder,"/",folder,"/",prefix,"_cleaned.nii",sep="")
nii_cleaned<-paste(folder,prefix,"_cleaned.nii",sep="")

data=data0[data0$V1 %in% dencl$V1,];
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
threshold<-(hl & dencl$V3>listminprob[cl])

#avtab<-read.table(paste(base_folder,"mip.dat",sep=""))
avtab<-read.table(paste(folder,prefix,"_mip.dat",sep=""))

dimvec<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
dimvec<-as.numeric(strsplit(dimvec," ")[[1]])

hc <- sub('FF$','77',heat.colors(nEP))
nx=dimvec[1]
ny=dimvec[2]
y=cenfile$V1+1
x=cenfile$V2+1

ord<-order(EPlist)

sizes<-rep(0,max(cl));
sizes[as.numeric(names(table(cl[threshold])))]=as.numeric(table(cl[threshold]));sizes_cl<-sizes;
ordcl<-order(sizes,decreasing=T);

#centers are the average response of each epoch for each cluster
centers<-matrix(nrow=max(cl),ncol=nEP);
centers_sd<-matrix(nrow=max(cl),ncol=nEP);
colors=rainbow(max(cl));

#create data frame of output.dat for cells in a given cluster above threshold, create mean and sd of cells in cluster
for(i in 1:max(cl)){
	tmpframe=data[cl==i & threshold,1:nEP+1];
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

resc <- function(x,zero=0.2){
	return(zero+(0.999-zero)*(x-min(x))/(max(x)-min(x)))
}

colnames(centers)<-EPlist
colnames(centers_sd)<-EPlist
write.table(x=data.frame(PREFIX=prefix,SIZE=sizes,centers),file=paste(folder,prefix,"_centers.dat",sep=""),row.names=F)
write.table(x=data.frame(PREFIX=prefix,SIZE=sizes,centers_sd),file=paste(folder,prefix,"_centers_sd.dat",sep=""),row.names=F)

plot_traj<-function(X,min=0,max=1,NS=1,AS=1){
	par(las=2,mar=c(3,10,2,3))
	ord2<-order(centers[X,]);
	pos<-barplot(centers[X,ord2],names.arg=EPlist[ord2],horiz=T,xlim=c(0,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS);
	arrows(y0=pos,y1=pos,x0=centers[X,ord2]-centers_sd[X,ord2],x1=centers[X,ord2]+centers_sd[X,ord2],code=3,length=0);
	abline(v=1);
}

plot_dp <- function(){
	plot(log(dencl$V4)[threshold],log(dencl$V3)[threshold],col=colors[ordcl[cl[threshold]]],pch=19);
	lines(rep(thsum[1],100),seq(-10,10,length=100),lty=2,col="green")
	lines(seq(-10,10,length=100),-nEP*seq(-10,10,length=100)+thsum[2],lty=2,col="green")
}

plot_dp_free <- function(c){
	logdist<-seq(min(log(dencl$V4)),max(log(dencl$V4)),length.out=100);
	logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
	colors<-rep("black",nrow(dencl));
	colors[cl==c]="red";
	plot(log(dencl$V4),log(dencl$V3),col=colors,pch=18);
	lines(rep(thsum[1],100),logprob,lty=2,col="green")
	lines(logdist,-nEP*logdist+thsum[2],lty=2,col="green")
}

plot_cluster<-function(X,P=19,S=1.5){
	colors=rainbow(max(cl));
	image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
	
	counter=1;
	for(i in X){
		points(x[cl==i & threshold],y[cl==i & threshold ],col=rgb(1,0,0,resc(dencl$V3[cl==i & threshold])),pch=P,cex=S)
		counter=counter+1;
	}
}

plot_cluster_h<-function(X,P=19,S=1){

#colors=c("#FF0000FF","#00FF00FF","#00FFFFFF","#0000FFFF","#FF00FFFF","#FFFF00FF")
	colors=rainbow(max(cl));
	image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=grey.colors(100),xaxt='n',yaxt='n',ann=FALSE)
	for(i in X) points(x[cl==(i) & threshold ],y[cl==(i) & threshold],col=colors[i],pch=P,cex=S)

}

plot_all<-function(X){
	image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=hc)
	for(i in X) points(x,y)
}

plot_ct <- function(c,ylim=c(0,1),maxlines=10^5){
	plot(rep(0,nEP),cex=0,ylim=ylim,xaxt='n',yaxt='n',xlab='',ylab="");
	counter=0;
	tmp<-data[cl==c,];
	tmpcl<-dencl[cl==c,];
	tmp<-tmp[order(tmpcl$V3,decreasing=T),];
	colalpha<-resc(sort(tmpcl$V3,decreasing=T))
	for(i in tmp$V1){
		if(counter<maxlines){
			lines(as.numeric(data[data$V1==i,1+1:nEP]),col=rgb(1,0,0,colalpha[counter+1]));
		}
		counter=counter+1;
	}
	lines(centers[c,],lwd=4);
}

plot_mix <- function(X,c=.4,NS,AS){
	layout(matrix(c(1,1,2,1,1,2,3,3,4),ncol=3));
	par(mar=c(3,4,1,1));
	plot_cluster(ordcl[X]);
	plot_ct(ordcl[X],c(-5,5));
	plot_traj(ordcl[X])
	par(mar=c(3,4,1,1));
	plot_dp_free(ordcl[X]);
}

A<-as.matrix(data[,1+1:nEP])
SVD<-svd(A);
plot_svd <- function(X,min=0,max=1,NS=1,AS=1){
	par(las=2,mar=c(3,20,2,3))
	ord2<-order(SVD$v[,X]);
	pos<-barplot(SVD$v[ord2,X],names.arg=EPlist[1:20][ord2],horiz=T,xlim=range(SVD$v[,X]),cex.axis=AS,cex.names=NS);
}


Fscale <- function(x){
	    (x-min(x))/(max(x)-min(x))
}




HM<-function(X,P=19,S=1){
	par(mfrow=c(2,2))
	a<-grep(paste("WIGGLY",X,"_",sep=""),EPlist)+1
		for(i in 1:length(a)){
		image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=grey.colors(100),xaxt='n',yaxt='n',ann=FALSE);
		title(main=EPlist[a[i]-1],cex.main=0.75)
		points(x,y,col=rgb(1,0,0,Fscale(data0[,a[i]])),pch=P,cex=S);
	}
}



