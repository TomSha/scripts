prefix="180607_PVN_barrage4_F1_1"

folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",prefix,sep="")
epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="")
nifti=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");


outputRCT<-read.table(paste(folder,"_output_rawCT.dat",sep=""))
output<-read.table(paste(folder,"_outputC.dat",sep=""))
noisethresh<-read.table(paste(folder,"_noise_thresh.dat",sep=""))
cenfile<-read.table(paste(folder,"_cell_centers.dat",sep=""))[noisethresh[,1],]
avtab<-read.table(paste(folder,"_mip.dat",sep=""))
dimvec<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
dimvec<-as.numeric(strsplit(dimvec," ")[[1]])
EPlength<-read.table(paste(folder,"_EPtimeC.dat",sep=""))$V1

nx=dimvec[1]
ny=dimvec[2]
y=cenfile$V1+1
x=cenfile$V2+1

ep<-read.table(epfile);
nEP<-nrow(ep)/2
ep$V1<-floor((ep$V1+90)*20);

epochs<-ep[,2][2:nEP*2]
#\\D not digit, [^BD] other than B & D
ord<-order(paste(gsub("\\D","",epochs),gsub("[^BD]","",epochs)))
epochs<-epochs[ord]
outputRCT[,2:ncol(outputRCT)]<-outputRCT[,2:ncol(outputRCT)][,ord]
EPlength<-EPlength[ord]

#minus 1 because we need to exclude the concentric
split<-(nEP-1)/2

#dot-bar/total
dotS_SF<-(outputRCT[,2:ncol(outputRCT)][,1:split*2]-outputRCT[,2:ncol(outputRCT)][,1:split*2-1])/(outputRCT[,2:ncol(outputRCT)][,1:split*2]+outputRCT[,2:ncol(outputRCT)][,1:split*2-1])
dotS90<-rowSums(dotS_SF[,1:5])/5
dotS270<-rowSums(dotS_SF[,6:10])/5
dotS<-rowSums(dotS_SF[,1:10])/10
dotselec<-data.frame(dotS_SF,dotS90,dotS270,dotS)

EPlist<-vector("list",length=4)
Olist<-vector("list",length=4)

EPlist[[1]]<-epochs[grep("BAR.270",x=epochs)]
Olist[[1]]<-outputRCT[,c(1,grep("BAR.270",x=epochs)+1)]
EPlist[[2]]<-epochs[grep("BAR.90",x=epochs)]
Olist[[2]]<-outputRCT[,c(1,grep("BAR.90",x=epochs)+1)]
EPlist[[3]]<-epochs[grep("DOT.270",x=epochs)]
Olist[[3]]<-outputRCT[,c(1,grep("DOT.270",x=epochs)+1)]
EPlist[[4]]<-epochs[grep("DOT.90",x=epochs)]
Olist[[4]]<-outputRCT[,c(1,grep("DOT.90",x=epochs)+1)]

EPlist<-lapply(EPlist,function(x) x[c(5,1,2,3,4)])
Olist<-lapply(Olist,function(x) x[,c(1,6,2,3,4,5)])

barDS_SF<-(Olist[[2]][,2:ncol(Olist[[2]])]-Olist[[1]][,2:ncol(Olist[[1]])])/(Olist[[2]][,2:ncol(Olist[[2]])]+Olist[[1]][,2:ncol(Olist[[1]])])
dotDS_SF<-(Olist[[4]][,2:ncol(Olist[[4]])]-Olist[[3]][,2:ncol(Olist[[3]])])/(Olist[[4]][,2:ncol(Olist[[4]])]+Olist[[3]][,2:ncol(Olist[[3]])])
barDS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("BAR_90",epochs)])-sum(x[grep("BAR_270",epochs)]))/sum(x[grep("BAR",epochs)]))
dotDS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("DOT_90",epochs)])-sum(x[grep("DOT_270",epochs)]))/sum(x[grep("DOT",epochs)]))
DS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("90",epochs)])-sum(x[grep("270",epochs)]))/sum(x))
DSdf<-data.frame(barDS_SF,dotDS_SF,barDS,dotDS,DS)


DS.mod<-lm(barDS ~ dotDS)

dotDSthresh<-dotDS>0.5|dotDS<(-0.5)
barDSthresh<-barDS>0.5|barDS<(-0.5)
dotS90thresh<-dotS90>0.5|dotS90<(-0.5)
dotS270thresh<-dotS270>0.5|dotS270<(-0.5)


globals<-data.frame(dotDS,barDS,dotS90,dotS270);
globalsthresh<-data.frame(dotDSthresh,barDSthresh,dotS90thresh,dotS270thresh);

#look at selectivity comparing different size dots. Compare 5 vs 17,26 and 30 degrees. Compare 10 vs 17,26,30 degrees
sizes<-as.integer(substr(EPlist[[3]],9,10))
sizeList<-vector("list",length=length(sizes))
fivedeg<-vector("list",length=3)
tendeg<-vector("list",length=3)
for (i in 1:length(sizeList)) sizeList[[i]]<-outputRCT[,2:ncol(outputRCT)][grepl("DOT",epochs)&grepl(sizes[i],epochs)]
for (i in 3:length(sizeList)) fivedeg[[i-2]]<-(sizeList[[1]]-sizeList[[i]])/(sizeList[[1]]+sizeList[[i]])
for (i in 3:length(sizeList)) tendeg[[i-2]]<-(sizeList[[2]]-sizeList[[i]])/(sizeList[[2]]+sizeList[[i]])

min_max<-function(x) (x-min(x))/(max(x)-min(x))
EPscale<-outputRCT[,2:ncol(outputRCT)]/rowMeans(outputRCT[,2:ncol(outputRCT)])
EPscale<-apply(EPscale,2,minMax)
EPscale270<-rowMeans(outputRCT[,2:ncol(outputRCT)][,grep("270",epochs)])
EPscale270<-minMax(EPscale270)
EPscale90<-rowMeans(outputRCT[,2:ncol(outputRCT)][,grep("90",epochs)])
EPscale90<-minMax(EPscale90)



#Functions--------------------------------------------------------------------------------------------------------------
#1 plot_cell(cell)		Plots the calcium trajectory of one cell
#2 plot_dotselec()		Plots a histogram of the selectivity of dots vs bars for all sizes and directions
#3 plot_global()		Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions 
#4 plot_corr()			Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
#5 plot_hist_resp()		Plots histogram of responses for all stimuli and a summary of the median response
#6 plot_sizeselec()		Plots a histogram of size selectivities
#7 plot_cellloc(cell)	Plots the location of a cell in the tectum
#8 plot_tuning(cell)	Plots a cell's response, dot vs bar selectivity, direction selectivity and location in the tectum

#1 Plots the calcium trajectory of one cell
plot_cell <- function(id,xmax=100000){

    par(mar=c(8,3,2,1))
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder,sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix,id),xaxt='n',ann=F);
	axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2)
	rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(xmax,nEP*2),col=rgb(0,1,0,.05))
}

#2 Plots the histogram of the selectivity of dots vs bars for all sizes and directions
plot_dotselec<-function(){
	par(mfrow=c(3,5));

	for (i in 1:ncol(dotselec)){
		hist(dotselec[,i],breaks=seq(min(dotselec),max(dotselec),l=15),xlim=c(-1,1),ylim=c(0,110),main=paste(epochs[i*2],"/",epochs[i*2-1],sep=""))
	}
}

#3 Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions
plot_global<-function(P=19,S=1.5,thresh=0){
	par(mfrow=c(2,4));
	hist(dotDS,xlim=c(-1,1),ylim=c(0,50),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="90 degrees = +1 & 270 degrees = -1");
	hist(barDS,xlim=c(-1,1),ylim=c(0,50),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="90 degrees = +1 & 270 degrees = -1");
	hist(dotS90,xlim=c(-1,1),ylim=c(0,35),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (90 degrees)");
	hist(dotS270,xlim=c(-1,1),ylim=c(0,35),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (270 degrees)");

	
	for(i in 1:4){
		globalsthresh[,i]<-globals[,i]>thresh|globals[,i]<(-thresh)
		colPal<-colorRampPalette(colors=c("blue","white","red"))
		cols<-colPal(50)[cut(globals[,i],breaks=50,labels=F)];
		image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
		points(x=cenfile[,2][globalsthresh[,i]],y=cenfile[,1][globalsthresh[,i]],col=cols,pch=P,cex=S)
	}
}

#4 Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
plot_corr<-function(){
	par(mfrow=c(1,2));
	plot(barDS ~ dotDS,xlim=c(-1,1),ylim=c(-1,1),main="bar DS vs dot DS")
	abline(DS.mod)
	plot(x=DS,y=dotS,xlim=c(-1,1),ylim=c(-1,1),main="direction selectivity vs dot selectivity")
}

#5 Plots histogram of responses for all stimuli and a summary of the median response
plot_hist_resp<-function(){
	#initialises dataframe
	outputnorm<-outputRCT

	for (i in 2:ncol(outputRCT)) outputnorm[,i]<-outputRCT[,i]/EPlength[i-1];
	par(mfrow=c(5,4));
	for (i in 2:ncol(outputRCT)){
		hist(outputnorm[,i],main=epochs[i-1],xlim=c(0,60),ylim=c(0,250),breaks=seq(min(outputnorm[,2:ncol(outputnorm)]),max(outputnorm[,2:ncol(outputnorm)]),l=30))
	}
	X11(); par(mar=c(8,3,2,1));

	cent<-barplot(apply(outputnorm[,2:ncol(outputnorm)],2,median),xaxt="n",main="response to each stimulus, normalised to epoch length");
	axis(side=1,labels=epochs,at=cent,las=2,cex.axis=0.75)
}

#6 Plots a histogram of size selectivities
plot_sizeselec<-function(){
	lmat=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T)
	layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5))
	minVal<-min(c(unlist(fivedeg),unlist(tendeg)))
	maxVal<-max(c(unlist(fivedeg),unlist(tendeg)))
	for (i in 1:3){
		for (j in 1:2) hist(fivedeg[[i]][,j],xlim=c(-1,1),ylim=c(0,70),breaks=seq(minVal,maxVal,l=20),main=paste(sizes[i+2]," degrees = +1 & 5 degrees = -1",sep=""))
	}
	
	X11(); layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5))
	for (i in 1:3){
		for (j in 1:2) hist(tendeg[[i]][,j],xlim=c(-1,1),ylim=c(0,70),breaks=seq(minVal,maxVal,l=20),main=paste(sizes[i+2]," degrees = +1 & 10 degrees = -1",sep=""))
	}
}


#7 Plots the location of a cell in the tectum
plot_cellloc<-function(cell){
	image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
#	points(x[outputRCT[cell,1]+1],y[outputRCT[cell,1]+1],col="red",pch=19,cex=1.5)
	points(x[cell],y[cell],col="red",pch=19,cex=1.5)
}

#8 Plots a cell's response, dot vs bar selectivity, direction selectivity and location in the tectum
plot_tuning<-function(cell){
	lmat=matrix(c(1,2,3,4,5,6,7,7),ncol=2,byrow=T)
	layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5))
	for (i in 1:length(EPlist)){

		center<-barplot(as.numeric(Olist[[i]][cell,2:ncol(Olist[[i]])]),ylim=c(0,max(outputRCT[cell,2:ncol(outputRCT)])),main=gsub("_5",EPlist[[i]][1],replacement=""))
		axis(1,at=center,labels=EPlist[[i]],cex.axis=1,las=2)

	}	
		label1=c(paste(epochs[1:split*2],"/",epochs[1:split*2-1],sep=""),"DOT_90/BAR_90","DOT_270/BAR_270","DOT/BAR")
		label2=c(paste(EPlist[[2]],"-",EPlist[[1]],sep=""),paste(EPlist[[4]],"-",EPlist[[3]],sep=""),"dotDS","barDS","DS")
		col1=c(paste(rep("green4",ncol(dotS_SF)/2)),rep("blue",ncol(dotS_SF)/2),rep("red",2),rep("black",1))
		col2=c(paste(rep("green4",ncol(barDS_SF))),rep("blue",ncol(dotDS_SF)),rep("red",2),rep("black",1))

		dotchart(as.numeric(dotselec[cell,]),labels=label1,col=col1,pch=19,xlim=c(-1,1),main="dot/bar selectivity")
		dotchart(as.numeric(DSdf[cell,]),labels=label2,col=col2,pch=19,xlim=c(-1,1),main="Direction selectivity")

		plot_cell(outputRCT[cell,1]+1)
		
		dc<-dev.cur()

		X11();
		plot_cellloc(cell);
		dev.set(dc)
}

plot_epoch_resp<-function(P=19,C=1.5){
	#image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
	#points(x,y,col=rgb(1,0,0,outputnorm[,epoch]),pch=P,cex=C)
	lmat=(matrix(1:20,ncol=4))
	layout(lmat)
	par(mar=c(0,0,2,0)); par(oma=c(0,0,2,0));
	Yflip<-rev(range(y))

	for (i in 1:length(epochs)){
		plot(x,y,col=rgb(0,0,0,EPscale[,i]),pch=P,cex=C,ann=T,xaxt="n",yaxt="n",bty="n",main=epochs[i],xlab="",ylab="",ylim=Yflip)
	}
	mtext("response heatmap of stimuli",outer=T)

	X11();
	lmat=matrix(1:2,ncol=2)
	layout(lmat)
	par(mar=c(0,0,2,0));par(oma=c(0,0,2,0));
	plot(x,y,col=rgb(0,0,0,EPscale270),pch=P,cex=C,ann=T,xaxt="n",yaxt="n",bty="n",main="270 degrees",xlab="",ylab="",ylim=Yflip,cex.main=0.9)
	plot(x,y,col=rgb(0,0,0,EPscale90),pch=P,cex=C,ann=T,xaxt="n",yaxt="n",bty="n",main="90 degrees",xlab="",ylab="",ylim=Yflip,cex.main=0.9)
	mtext("response heatmap of stimuli moving in 270 or 90 degrees direction\n",outer=T)

}


