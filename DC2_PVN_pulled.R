#martincolorscale=c("#6ef914ff","#060cf1ff","#f106e3ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");
martincolorscale=c("green4","royalblue4","#f3fa54ff","#1f1c25ff","deeppink4","#176462ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff","#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff");
#args <- commandArgs(trailingOnly = TRUE)

#prefix<- args[1]

#CHANGED 1 removed min density threshold for halo removal

require(plotrix)
require(mmand)
prefix_list<-read.table("/home/meyer-lab/ToolBox_v3.6/scripts/prefix_list_bars_dots")$V1	                					#CHANGE HERE
results_toolbox<-"results_toolbox_sizes"																					#ALSO CHANGE HERE
folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/pulled_results_bars_dots",sep="") 	#ALSO CHANGE HERE	
niftifolder="/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/"
epfile=paste(folder,"/pulled_epochsC.dat",sep="");														
#how many nearest neighbours need to be in same cluster for cell not to be removed as a halo
KNN=3
EP<-read.table(epfile);
nEP=nrow(EP)
EPlist<-vector(mode="character", length=nEP);
EPlist<-as.character(EP$V1[1:nEP])


thsum<-as.numeric(read.table(paste(folder,"/pulled_thresholds.dat",sep="")));
data=read.table(paste(folder,"/all_centers.dat",sep=""));
dencl<-read.table(paste(folder,"/pulled_summary.dat",sep=""));
KNNR<-read.table(paste(folder,"/pulled_KNNrecords.dat",sep=""));

inc_list<-vector("list",length(prefix_list));
noisethresh_list<-vector("list",length(prefix_list));
dimvec_list<-vector("list",length(prefix_list))
center_list<-vector("list",length(prefix_list));
center_all_list<-vector("list",length(prefix_list));
avtab_list<-vector("list",length(prefix_list));
dencl_list<-vector("list",length(prefix_list));
corvec_list<-vector("list",length(prefix_list));
ep_list<-vector("list",length(prefix_list));
epochsC_list<-vector("list",length(prefix_list));
output_list<-vector("list",length(prefix_list));
output_raw_list<-vector("list",length(prefix_list));

for(i in 1:length(prefix_list)){
	prefix=prefix_list[i]
	nifti=paste(niftifolder,"/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");
	
	dencl_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_summary.dat",sep=""))#CHANGED HERE!!!!!
	dencl_list[[i]]<-dencl_list[[i]][order(dencl_list[[i]]$V1),];
	inc_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_inc.dat",sep=""));
	noisethresh_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_noise_thresh.dat",sep=""))
#	center_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_cell_centers.dat",sep=""))[noisethresh_list[[i]]$V1,][dencl_list[[i]]$V1+1,][inc_list[[i]][,1],];
	center_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_cell_centers.dat",sep=""))#[noisethresh_list[[i]]$V1,][inc_list[[i]][,1],];
	center_all_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_cell_centers.dat",sep=""))
	dencl_list[[i]]<-dencl_list[[i]][inc_list[[i]][,1],];
	dimvec_list[[i]]<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dimvec_list[[i]]<-as.numeric(strsplit(dimvec_list[[i]]," ")[[1]])
	avtab_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_mip.dat",sep=""))
	corvec_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_corvec.dat",sep=""))
	ep_list[[i]]<-read.table(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep=""))
	ep_list[[i]]$V1<-(ep_list[[i]]$V1+90)*20
	epochsC_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_epochsC.dat",sep=""))$V1
#	output_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_outputCTM.dat",sep=""))
	ord<-paste(gsub("\\D","",epochsC_list[[i]]),gsub("[^BD]","",epochsC_list[[i]]))#new
	ord<-order(gsub("5","05",ord))#new NB is overwritten later on
#	output_list[[i]][,2:ncol(output_list[[i]])]<-output_list[[i]][,2:ncol(output_list[[i]])][,ord]
	output_raw_list[[i]]<-read.table(paste(niftifolder,"/",prefix,"/im/im_slice_1/",results_toolbox,"/",prefix,"_output_rawCT.dat",sep=""))
	ep<-ep_list[[i]]$V2[2:(nrow(ep_list[[i]])/2)*2]
	ep_ord<-paste(gsub("\\D","",ep),gsub("[^BD]","",ep))#new
	ep_ord<-order(gsub("5","05",ep_ord))#new
	ep<-ep[ep_ord]#new
#	output_raw_list[[i]][,2:ncol(output_raw_list[[i]])]<-output_raw_list[[i]][,2:ncol(output_raw_list[[i]])][,order(paste(gsub("\\D","",ep),gsub("[^BD]","",ep)))]
	output_raw_list[[i]][,2:ncol(output_raw_list[[i]])]<-output_raw_list[[i]][,2:ncol(output_raw_list[[i]])][,ep_ord]
}

total_nEP<-nrow(ep_list[[1]])/2
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
#threshold<-(hl & dencl$V3>listminprob[cl])     #CHANGED 1
threshold<-hl
hc <- sub('FF$','77',heat.colors(nEP))
write.table(file=paste(folder,"/pulled_inc.dat",sep=""),threshold)

ord<-paste(gsub("\\D","",EPlist),gsub("[^BD]","",EPlist))
ord<-order(gsub("5","05",ord))

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

dotS_list<-lapply(output_raw_list,function(output) (output[,2:ncol(output)][,1:(length(ep)/2)*2] - output[,2:ncol(output)][,1:(length(ep)/2)*2-1]) / (output[,2:ncol(output)][,1:(length(ep)/2)*2] + output[,2:ncol(output)][,1:(length(ep)/2)*2-1]))

#dotS_list<-lapply(output_raw_list, function(output_raw) (output_raw[,grep("DOT_270_30",ep)+1]-output_raw[,grep("BAR_270_30",ep)+1])/(output_raw[,grep("DOT_270_30",ep)+1]+output_raw[,grep("BAR_270_30",ep)+1]))

dotDS_list<-lapply(output_raw_list, function(output_raw) (rowSums(output_raw[,grep("DOT_90",ep)+1])-rowSums(output_raw[,grep("DOT_270",ep)+1]))/(rowSums(output_raw[,grep("DOT_90",ep)+1])+rowSums(output_raw[,grep("DOT_270",ep)+1])))

#dotDS_list<-lapply(output_raw_list, function(output_raw) (output_raw[,grep("DOT_90_5",ep)+1]-output_raw[,grep("DOT_270_5",ep)+1])/(output_raw[,grep("DOT_90_5",ep)+1]+output_raw[,grep("DOT_270_5",ep)+1]))


dot_ent_list<-lapply(output_raw_list, function(output) output[,2:ncol(output)]/rowSums(output[,2:ncol(output)]))
dot_ent_list<-lapply(dot_ent_list, function(ent) apply(ent,1,function(x) sum(x*log2(1/x))))
max_ent<-log2(length(ep))

#FUNCTIONS

plot_traj<-function(X,min=0,max=1,NS=1,AS=1,CC=martincolorscale[X+21]){
	par(las=2,mar=c(3,10,2,3))
	ord2<-order(centers[X,]);
	pos<-barplot(centers[X,ord2],names.arg=EPlist[ord2],horiz=T,xlim=c(0,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS,col=CC);
#	pos<-barplot(centers[X,ord2],names.arg=c("POSTERIOR","ANTERIOR"),horiz=T,xlim=c(-1,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS,col=CC);
	arrows(y0=pos,y1=pos,x0=centers[X,ord2]-centers_sd[X,ord2],x1=centers[X,ord2]+centers_sd[X,ord2],code=3,length=0,col="black");
}

plot_box<-function(X){
	par(mar=c(3,7,2,2))
#	ord2<-order(centers[X,]);
	fish<-sum(!table(data$V2[cl==X & threshold])==0)
#	boxplot(data[cl==X & threshold,4+ord],names=EPlist[ord],las=2,horizontal=T,main=paste("cluster",X,"in",fish,"/",length(prefix_list),"fish"),col=c(rep(martincolorscale[1],5),rep(martincolorscale[2],5)))
	boxplot(data[cl==X & threshold,4+ord],names=rep("",2),las=2,horizontal=T,main="",col=c(rep(martincolorscale[4],5),rep(martincolorscale[2],5)))
#	abline(v=0, lty=2,col="blue")
}


plot_box_experiment<-function(cluster,experiment){
	par(mar=c(3,7,2,2))
	output<-output_list[[experiment]][,2:ncol(output_list[[experiment]])][(output_list[[experiment]][,1] %in% cluster_cells_ID(cluster)[[experiment]]),]
	if(nrow(output)==0){
		plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
	}else{
		boxplot(output,horizontal=T,names=EPlist[ord],las=2,col=c(rep(martincolorscale[1],5),rep(martincolorscale[2],5)))
	}
}

plot_dp <- function(){
	plot(log(dencl$V4)[threshold],log(dencl$V3)[threshold],col=colors[ordcl[cl[threshold]]],pch=19,cex=0.6,xlab="log(distance)",ylab="log(density)");
	logdist<-seq(min(log(dencl$V4)),range(dencl$V4,finite=T)[2],length.out=100);
	logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
#	lines(rep(thsum[1],100),logprob,lty=2,col="green")
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

plot_dp_cluster <- function(c){
	logdist<-seq(min(log(dencl$V4)),range(dencl$V4,finite=T)[2],length.out=100);
	logprob<-seq(min(log(dencl$V3)),max(log(dencl$V3)),length.out=100);
	colors<-rep("black",nrow(dencl));
	colors[cl[threshold]==c]=martincolorscale[11];
	plot(log(dencl$V4)[threshold],log(dencl$V3)[threshold],col=colors,pch=19,cex=1,xlab="log(distance)",ylab="log(density)",cex.lab=1.5);
	#lines(rep(thsum[1],100),logprob,lty=2,col="green")
	lines(logdist,-nEP*logdist+thsum[2],lty=2,col=martincolorscale[1])
}


plot_corvec<-function(P=19,S=1.5,DR=50){
	lmat=matrix(1:8,ncol=2)
	layout(lmat)
	par(mar=c(0.5,0.5,0.5,0.5))
	colPal<-colorRampPalette(colors=c("blue","white","red"))
	

	for(k in 1:length(prefix_list)){
		nx=dimvec_list[[k]][1]
		ny=dimvec_list[[k]][2]
		image(1:nx,1:ny,matrix(avtab_list[[k]]$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE);
		cols<-colPal(50)[cut(corvec_list[[k]][,1],breaks=50,labels=F)];
		points(y=center_all_list[[k]][,1],x=center_all_list[[k]][,2],col=cols,pch=P,cex=S);
	}
}

plot_image<-function(cluster,experiment,P=19,S=1.5,thresh=T,cell="ALL",c=1){
#dencl_list[[k]]$V6 : C1 ID from fish k
#data$V4			: C1 ID
#data$V2			: prefix

	y=center_list[[experiment]]$V1+1
	x=center_list[[experiment]]$V2+1
	nx=dimvec_list[[experiment]][1]
	ny=dimvec_list[[experiment]][2]
	image(1:nx,1:ny,matrix(avtab_list[[experiment]]$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE);
	for (i in 1:length(cluster)){
		if(thresh){
		selec=dencl_list[[experiment]]$V6 %in% (data$V4[data$V2==prefix_list[experiment] & cl==cluster[i] & threshold]-1)
		}else{
		selec=dencl_list[[experiment]]$V6 %in% (data$V4[data$V2==prefix_list[experiment] & cl==cluster[i]]-1)
		}
		if(cell=="ALL"){
		points(x[ selec  ],y[selec ],col=martincolorscale[c],pch=P,cex=S)
		}else{
		points(x[ selec  ][cell],y[selec ][cell],col=martincolorscale[i],pch=P,cex=S)
		}	
	}
}



cluster_cells_ID<-function(cluster){
	mapply(function(dencl_list,prefix) dencl_list$V1[dencl_list$V6 %in% (data$V4[data$V2==prefix & cl==cluster & threshold]-1)], dencl_list,prefix_list)
}

read_trace<-function(experiment,id){	
	ep<-ep_list[[experiment]]
	prefix<-prefix_list[experiment]
	bin<-paste(niftifolder,prefix,"/im/im_slice_1/results_toolbox/",prefix,sep="")
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",(id+1)," ",bin,sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	returnList<-list("prefix"=prefix,"ep"=ep,"trace"=t_dt)
}

plot_trace<-function(cluster,experiment,cell,xmax=10000){
	id<-cluster_cells_ID(cluster)[[experiment]][cell]
	params<-read_trace(experiment,id)
	par(mar=c(8,3,2,1))
	plot((params$trace),type='l',ylim=(c(0,max(params$trace))),main=paste(params$prefix,id),xaxt='n',ylab="",xlab="");
	axis(1,at=params$ep$V1[1:total_nEP*2-1],label=params$ep$V2[1:total_nEP*2-1],las=2);
	rect(xleft=params$ep$V1[1:total_nEP *2 -1],xright=params$ep$V1[1:total_nEP *2 ],ybottom=rep(-1000,total_nEP*2),ytop=rep(xmax,total_nEP*2),col=rgb(0,1,0,.05));
	returnList<-list("cluster"=cluster,"prefix"=prefix,"cell_id"=id)
	return(returnList)
}

plot_trace_average<-function(cluster,experiment,xmax=10000){
	par(mar=c(8,3,2,2))
	cells<-cluster_cells_ID(cluster)[[experiment]]
	if(length(cells)==0){
	#if there are no cells in the cluster for this fish, plot empty plot	
		plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
	}else{
		params<-lapply(cells,function(cells) read_trace(experiment,cells))
		bin_cluster<-lapply(params,"[[",3)
		bin_cluster<-lapply(bin_cluster,function(x) x/mean(x))
		bin_cluster<-sapply(bin_cluster,unlist)
		plot((rowMeans(bin_cluster)),type='l',ylim=(c(0,max(rowMeans(bin_cluster)))),main=paste(cluster),xaxt='n',ylab="",xlab="",col=martincolorscale[2]);
		axis(1,at=params[[1]]$ep$V1[1:total_nEP*2-1],label=params[[1]]$ep$V2[1:total_nEP*2-1],las=2);
		rect(xleft=params[[1]]$ep$V1[1:total_nEP *2 -1],xright=params[[1]]$ep$V1[1:total_nEP *2 ],ybottom=rep(-1000,total_nEP*2),ytop=rep(xmax,total_nEP*2),col=rgb(0,0,0,.05));
	}
}


plot_cell<-function(cluster,experiment,cell){
	output<-output_list[[experiment]]
	lmat<-matrix(c(1,1,1,2,3,4),ncol=3,byrow=T)
	par(mar=c(5,3,2,1))
	layout(lmat)
	
	id<-plot_trace(cluster,experiment,cell)[3]
	barplot(as.numeric(output[output$V1==id,2:ncol(output)]),names.arg=EPlist[ord],las=2)
	plot_box(cluster)
	plot_image(cluster,experiment,cell=cell)	
}
	
plot_cluster<-function(cluster,new_plot=F){
	if (dev.cur()==1|new_plot){
		dev.new(width=18,height=12)
	}
	lmat=matrix(c(1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8,9,9,10,10,9,9,10,10),ncol=4,byrow=T);
	layout(lmat)
	colors=rainbow(max(cl));

	for(experiment in 1:length(prefix_list)){
		plot_trace_average(cluster,experiment)
		par(mar=c(1,1,1,1))
		plot_image(cluster,experiment)
	}
	plot_box(cluster)
	plot_dp_cluster(cluster);	

	dev<-dev.cur()

	if(length(dev.list())<2|new_plot){
		X11();
		par(mfrow=c(2,2)); par(oma=c(0,0,2,0));
		for (i in 1:length(prefix_list)) plot_box_experiment(cluster,i)
		mtext(paste("cluster",cluster),outer=T)
		dev.set(dev)
	}else{
		if(dev.cur()==2){
			dev.set(3)
		}else{
			dev.set(2)
		}
		par(mfrow=c(2,2))
		for (i in 1:length(prefix_list)) plot_box_experiment(cluster,i)
		dev.set(dev)
	}
}

plot_dotS_cluster<-function(cluster){
	index<-mapply(function(output,cells) output[,1] %in% cells, output_raw_list,cluster_cells_ID(cluster))
	index<-unlist(index)
#	dotS<-unlist(dotS_list)
	dotS<-do.call(rbind,dotS_list)[index,]
	dotS<-rowSums(dotS)/10
	hist(dotS,xlim=c(-1,1),breaks=seq(min(dotS),max(dotS),l=7),main="",col=martincolorscale[1],xlab="")
	legend("topleft",inset=c(0.05,0.01),legend=c("1 = dot selective","-1 = bar selective"),xpd=TRUE,cex=1)


}

plot_dotS_cluster_box<-function(cluster){
  index<-mapply(function(output,cells) output[,1] %in% cells, output_raw_list,cluster_cells_ID(cluster))
  index<-unlist(index)
  #	dotS<-unlist(dotS_list)
  dotS_cluster<-do.call(rbind,dotS_list)[index,]
  dotS_cluster<-rowSums(dotS_cluster)/10
  dotS<-do.call(rbind,dotS_list)
  dotS<-rowSums(dotS)/10
  boxplot(dotS_cluster,ylim=c(-1,1),main="",col=martincolorscale[1],xlab="",at=-0.5,wex=0.5)
  boxplot(dotS,add=T,at=0.5,wex=0.5)
  #legend("topleft",inset=c(0.05,0.01),legend=c("1 = dot selective","-1 = bar selective"),xpd=TRUE,cex=1)
  
  
}

plot_dotDS_cluster<-function(cluster){
	index<-mapply(function(output,cells) output[,1] %in% cells, output_raw_list,cluster_cells_ID(cluster))
	index<-unlist(index)
	DS<-as.numeric(unlist(dotDS_list)[index])
	hist(DS,xlim=c(-1,1),breaks=seq(min(DS),max(DS),l=7),main="",col=martincolorscale[1],xlab="")
	legend("topright",inset=c(0.05,0.01),legend=c("1 = 90ᵒ","-1 = 270ᵒ selective"),xpd=TRUE,cex=1)

}

plot_ent_cluster<-function(cluster){
	index<-mapply(function(output,cells) output[,1] %in% cells, output_raw_list,cluster_cells_ID(cluster))
	index<-unlist(index)
	ent<-as.numeric(unlist(dot_ent_list)[index])
	boxplot(x=ent,las=2,cex.axis=0.8,boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),ylim=c(0,max_ent),ylab="entropy")
	boxplot(ent,boxwex=0.5,add=T,at=0.75,yaxt="n",xaxt="n",boxfill=martincolorscale[1])
	boxplot(unlist(dot_ent_list),boxwex=0.5,add=T,at=1.25,yaxt="n",xaxt="n",boxfill=martincolorscale[2])
	abline(h=max_ent,lty=3)
	legend("bottomleft",inset=c(0.05,0.01),legend=c(paste("cluster",cluster),"all cells"),fill=martincolorscale[1:2],xpd=TRUE,cex=1.5)

}



plot_analysis<-function(cluster){
	dev.new(width=18,height=5)
	par(mfrow=c(1,4))
	plot_box(cluster)
	plot_ent_cluster(cluster)
	plot_dotS_cluster(cluster)
	plot_dotDS_cluster(cluster)
}

#plot_cluster<-function(cluster,P=19,S=1.5,nothresh=F,traj_device=NA){

#	if(nothresh){
#	for(i in 1:length(X)){
#		lmat=matrix(c(1,1,3,1,1,4,2,2,5,2,2,6),ncol=3,byrow=T);
#		layout(lmat)
#		colors=rainbow(max(cl));
#		plot_box(X[i])
#		plot_dp_cluster(X[i]);	
#		par(mar=c(0.5,0.5,0.5,0.5))
#		for(experiment in 1:length(prefix_list)){
#	    	plot_image(cluster,experiment)			
#			}
#		}
#		}else{
#	for(i in 1:length(X)){
#	lmat=matrix(c(1,1,3,1,1,4,2,2,5,2,2,6),ncol=3,byrow=T);
#		layout(lmat)
#		colors=rainbow(max(cl));
#		plot_box(X[i])
#		par(mar=c(5,7,2,2))
#		plot_dp_cluster(X[i]);	
#		par(mar=c(0.5,0.5,0.5,0.5))
#		for(k in 1:length(prefix_list)){
#	   		plot_image(cluster,experiment)
#		}
#		}
#		}
#	}


