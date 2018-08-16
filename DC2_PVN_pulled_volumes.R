#args <- commandArgs(trailingOnly = TRUE)

#prefix<- args[1]

require(plotrix)
require(mmand)
martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");
prefix_list<-as.character(read.table("/home/meyer-lab/ToolBox_v3.6/scripts/prefix_list_vol")[,1]) 
folder="/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/"

##CHANGE HERE!!!!####################
results="/results_toolbox/"
pulled="pulled_results_volumes"
#####################################

epfile=paste(folder,pulled,"/pulled_epochsC.dat",sep=""); #TO ADD in the folder
#how many nearest neighbours need to be in same cluster for cell not to be removed as a halo
KNN=3
EP<-read.table(epfile);
nEP=nrow(EP)
EPlist<-vector(mode="character", length=nEP);
EPlist<-as.character(EP$V1[1:nEP])

thsum<-as.numeric(read.table(paste(folder,pulled,"/pulled_thresholds.dat",sep=""))); 
data=read.table(paste(folder,pulled,"/all_centers.dat",sep=""));
dencl<-read.table(paste(folder,pulled,"/pulled_summary.dat",sep=""));
KNNR<-read.table(paste(folder,pulled,"/pulled_KNNrecords.dat",sep=""));

inc_list<-vector("list",length(prefix_list));
dimvec_list<-vector("list",length(prefix_list))
center_list<-vector("list",length(prefix_list));
dencl_list<-vector("list",length(prefix_list));
for(i in 1:length(prefix_list)){
	prefix=prefix_list[i]
	nifti=paste(folder,prefix,"/nii/rslice10/rim10.nii",sep="");

	dencl_list[[i]]<-read.table(paste(folder,prefix,results,prefix,"_summary.dat",sep=""))
	dencl_list[[i]]<-dencl_list[[i]][order(dencl_list[[i]]$V1),];
	inc_list[[i]]<-read.table(paste(folder,prefix,results,prefix,"_inc.dat",sep=""));
	center_list[[i]]<-read.table(paste(folder,prefix,results,prefix,"_cell_centers.dat",sep=""))[dencl_list[[i]]$V1+1,][inc_list[[i]][,1],];
	dencl_list[[i]]<-dencl_list[[i]][inc_list[[i]][,1],];
	dimvec_list[[i]]<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dimvec_list[[i]]<-as.numeric(strsplit(dimvec_list[[i]]," ")[[1]])
}

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
#threshold<-(hl & dencl$V3>listminprob[cl])
# Tom's idea is to release the second constraint on the minimum probability so 
threshold<-hl
hc <- sub('FF$','77',heat.colors(nEP))
write.table(file=paste(folder,pulled,"/pulled_inc.dat",sep=""),threshold)

ord<-order(EPlist)

sizes<-rep(0,max(cl));
sizes[as.numeric(names(table(cl[threshold])))]=as.numeric(table(cl[threshold]));
ordcl<-order(sizes,decreasing=T);

centers<-matrix(nrow=max(cl),ncol=nEP);
centers_sd<-matrix(nrow=max(cl),ncol=nEP);
colors=rainbow(max(cl));

for(i in 1:max(cl)){
	tmpframe=data[cl==i & threshold,1:nEP+4];
	tmpdencl=dencl[cl==i & threshold,];
	if(!is.null(nrow(tmpframe))){
		if(nrow(tmpframe)==0){
			centers[i,] <- rep(0,nEP);
			centers_sd[i,] <- rep(0,nEP);
		}
		if(nrow(tmpframe)>1){
			#centers[i,] <- apply(tmpframe,2,mean);
			centers[i,] <- as.numeric(tmpframe[which.max(tmpdencl$V3),])
			centers_sd[i,] <- apply(tmpframe,2,sd);
		}
		if(nrow(tmpframe)==1){
			centers[i,]<- as.numeric(tmpframe[1,]);
			centers_sd[i,]<-rep(0,nEP)
		}

	} else {
		centers[i,] <- tmpframe;
		centers_sd[i,] <- rep(0,nEP);
	}
}

plot_traj<-function(X,min=0,max=1,NS=1,AS=1,CC="blue"){
	par(las=2,mar=c(3,10,2,3))
	ord2<-order(centers[X,]);
#	pos<-barplot(centers[X,ord2],names.arg=EPlist[ord2],horiz=T,xlim=c(0,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS,col=CC);
	pos<-barplot(centers[X,ord2],names.arg=EPlist[ord2],horiz=T,xlim=c(-1,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS,col=CC,main=X[1]);
	arrows(y0=pos,y1=pos,x0=centers[X,ord2]-centers_sd[X,ord2],x1=centers[X,ord2]+centers_sd[X,ord2],code=3,length=0,col="black");
}

plot_dp <- function(){
	logdist<-seq(min(log(dencl$V4)),range(dencl$V4,finite=T)[2],length.out=100);
    logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
	plot(log(dencl$V4)[threshold],log(dencl$V3)[threshold],col=colors[ordcl[cl[threshold]]]);
	lines(rep(thsum[1],100),logprob,lty=2,col="green")
	lines(logdist,-nEP*logdist+thsum[2],lty=2,col="green")
}

plot_dp_free <- function(c){
	logdist<-seq(min(log(dencl$V4)),range(dencl$V4,finite=T)[2],length.out=100);
	logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
	colors<-rep("black",nrow(dencl));
	colors[cl==c]="red";
	plot(log(dencl$V4),log(dencl$V3),col=colors,pch=18);
	lines(rep(thsum[1],100),logprob,lty=2,col="green")
	lines(logdist,-nEP*logdist+thsum[2],lty=2,col="green")
}

plot_cluster<-function(X,P=19,S=1.5,nothresh=F,traj_device=NA){
	counter=0;
	for(k in 1:length(prefix_list)){
	    y=center_list[[k]]$V1
		x=center_list[[k]]$V2
		z=center_list[[k]]$V3

		nx=dimvec_list[[k]][1]
		ny=dimvec_list[[k]][2]
		dftmp<-list()
		for(CLID in 1:max(cl)){
			selec=dencl_list[[k]]$V6 %in% (data$V4[data$V2==prefix_list[k] & cl==CLID & threshold]-1)
#			selec=dencl_list[[k]]$V6 %in% (data$V4[data$V2==prefix_list[k] & cl==CLID ]-1)
			if(sum(selec)>0) dftmp[[CLID]] <- data.frame(x[selec],y[selec],z[selec],CLID, 1)
		}

		dftmpbind=do.call(rbind,dftmp)
		pointfile=paste("points3D_",prefix_list[k],".dat",sep="")
		write.table(dftmpbind,pointfile,row.names=F,col.names=F)
		system(paste("../bin/3Dplot", prefix_list[k],pointfile,500*(counter%%3),500*floor(counter/3),paste(X,collapse=" "),"&"))
		counter=counter+1
	}
}

plot_cluster_fun<-function(custom_cl,P=19,S=1.5,nothresh=F,traj_device=NA){
	lmat=matrix(1:8,ncol=2);
	layout(lmat)
	par(mar=c(.5,.5,.5,.5))
	colors=rainbow(max(cl));
	counter=0;
	for(k in 1:length(prefix_list)){
	    y=center_list[[k]]$V1
		x=center_list[[k]]$V2
		z=center_list[[k]]$V3

		nx=dimvec_list[[k]][1]
		ny=dimvec_list[[k]][2]
		dftmp<-list()

		for(CLID in 0:1){
			selec=dencl_list[[k]]$V6 %in% (data$V4[data$V2==prefix_list[k] & custom_cl==CLID & threshold]-1)
			if(sum(selec)>0) dftmp[[CLID+1]] <- data.frame(x[selec],y[selec],z[selec],CLID)
		}

		dftmpbind=do.call(rbind,dftmp)
		write.table(dftmpbind,file=paste("points3D",prefix_list[k],".dat",sep=""),row.names=F,col.names=F)
		system(paste("../bin/3Dplot",paste("points3D",prefix_list[k],".dat",sep=""),500*(counter%%3),500*floor(counter/3),1,"&"))
		counter=counter+1
	}
}


