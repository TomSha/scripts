prefix="180727_PVN_barrage4_v10_F1_1"

folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox/",prefix,sep="")
epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="")
niftifolder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/",sep="");

EPlength<-read.table(paste(folder,"_rslice10_EPtimeC.dat",sep=""))$V1
outputRCT<-read.table(paste(folder,"_combined_output_rawCT.dat",sep=""))
output<-read.table(paste(folder,"_combined_outputCT.dat",sep=""))
noisethresh<-read.table(paste(folder,"_noise_thresh.dat",sep=""))
cenfile<-read.table(paste(folder,"_cell_centers.dat",sep=""))[noisethresh[,1],]
#avtab<-read.table(paste(folder,"_mip.dat",sep=""))
#dimvec<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
#dimvec<-as.numeric(strsplit(dimvec," ")[[1]])

#nx=dimvec[1]
#ny=dimvec[2]
#y=cenfile$V1+1
#x=cenfile$V2+1

minSlice<-min(outputRCT[,1])
numSlices<-max(outputRCT[,1])-minSlice
outputS_list<-vector("list",length=numSlices);
for (i in 0:numSlices){
	slice<-i+minSlice
	outputS_list[[i+1]]<-read.table(paste(folder,"_rslice",slice,"_output_rawCT.dat",sep=""))
	outputS_list[[i+1]][,2:ncol(outputS_list[[i+1]])]<-outputS_list[[i+1]][,2:ncol(outputS_list[[i+1]])][,ord]
	}


ep<-read.table(epfile);
nEP<-(nrow(ep)/2)-1
ep$V1<-floor((ep$V1+90)*4);

epochs<-ep[3:nrow(ep),2][1:nEP*2]
#\\D not digit, [^BD] other than B & D
ord<-order(paste(gsub("\\D","",epochs),gsub("[^BD]","",epochs)))
epochs<-epochs[ord]
outputRCT[,3:ncol(outputRCT)]<-outputRCT[,3:ncol(outputRCT)][,ord]
EPlength<-EPlength[ord]

#minus 1 because we need to exclude the concentric
split<-(nEP)/2

#dot-bar/total
dotS_SF<-(outputRCT[,3:ncol(outputRCT)][,1:split*2]-outputRCT[,3:ncol(outputRCT)][,1:split*2-1])/(outputRCT[,3:ncol(outputRCT)][,1:split*2]+outputRCT[,3:ncol(outputRCT)][,1:split*2-1])
dotS90<-rowSums(dotS_SF[,1:5])/5
dotS270<-rowSums(dotS_SF[,6:10])/5
dotS<-rowSums(dotS_SF[,1:10])/10
dotselec<-data.frame(dotS_SF,dotS90,dotS270,dotS)


#calculate direction selectivity (for bars and dots)
EPsplit<-c("BAR.270","BAR.90","DOT.270","DOT.90")
reord<-c(5,1,2,3,4)

EPlist<-vector("list",length=length(EPsplit))
Olist<-vector("list",length=length(EPsplit))
for (i in 1:length(EPsplit)){
    EPlist[[i]]<-epochs[grep(EPsplit[i],x=epochs)]
    Olist[[i]]<-outputRCT[,c(1,2,grep(EPsplit[i],x=epochs)+2)]
}
EPlist<-lapply(EPlist,function(x) x[reord])
Olist<-lapply(Olist,function(x) x[,c(1,2,(reord+2))])


barDS_SF<-(Olist[[2]][,3:ncol(Olist[[2]])]-Olist[[1]][,3:ncol(Olist[[1]])])/(Olist[[2]][,3:ncol(Olist[[2]])]+Olist[[1]][,3:ncol(Olist[[1]])])
dotDS_SF<-(Olist[[4]][,3:ncol(Olist[[4]])]-Olist[[3]][,3:ncol(Olist[[3]])])/(Olist[[4]][,3:ncol(Olist[[4]])]+Olist[[3]][,3:ncol(Olist[[3]])])
barDS<-apply(outputRCT[,3:ncol(outputRCT)],1,function(x) (sum(x[grep("BAR_90",epochs)])-sum(x[grep("BAR_270",epochs)]))/sum(x[grep("BAR",epochs)]))
dotDS<-apply(outputRCT[,3:ncol(outputRCT)],1,function(x) (sum(x[grep("DOT_90",epochs)])-sum(x[grep("DOT_270",epochs)]))/sum(x[grep("DOT",epochs)]))
DS<-apply(outputRCT[,3:ncol(outputRCT)],1,function(x) (sum(x[grep("90",epochs)])-sum(x[grep("270",epochs)]))/sum(x))
DSdf<-data.frame(barDS_SF,dotDS_SF,barDS,dotDS,DS)


DS.mod<-lm(barDS ~ dotDS)
globals<-data.frame(dotDS,barDS,dotS90,dotS270);

#look at selectivity comparing different size dots. Compare 5 vs 17,26 and 30 degrees. Compare 10 vs 17,26,30 degrees
sizes<-as.integer(substr(EPlist[[3]],9,10))
sizeList<-vector("list",length=length(sizes))
fivedeg<-vector("list",length=3)
tendeg<-vector("list",length=3)
for (i in 1:length(sizeList)) sizeList[[i]]<-outputRCT[,3:ncol(outputRCT)][grepl("DOT",epochs)&grepl(sizes[i],epochs)]
for (i in 3:length(sizeList)) fivedeg[[i-2]]<-(sizeList[[1]]-sizeList[[i]])/(sizeList[[1]]+sizeList[[i]])
for (i in 3:length(sizeList)) tendeg[[i-2]]<-(sizeList[[2]]-sizeList[[i]])/(sizeList[[2]]+sizeList[[i]])


#look at the mean responses across epochs
#normalise across cells to mean and then across epochs to epoch length
cellNorm<-t(apply(outputRCT[,3:ncol(outputRCT)],1,function(x) x/mean(x)))
for (j in 1:ncol(cellNorm)){
        cellNorm[,j]<-cellNorm[,j]/EPlength[j]
 }


EPMeans<-colMeans(cellNorm)
#applies sd over epochs for each fish
EPSD<-apply(cellNorm,2,sd)

EPMean_list<-vector("list",length=length(EPsplit))
EPSD_list<-vector("list",length=length(EPsplit))
CoVa_list<-vector("list",length=length(EPsplit))
for (i in 1:length(EPsplit)){
	    EPMean_list[[i]]<-EPMeans[grep(EPsplit[i],x=epochs)]
		    EPSD_list[[i]]<-EPSD[grep(EPsplit[i],x=epochs)]
			    CoVa_list[[i]]<-EPSD_list[[i]]/EPMean_list[[i]]
}

EPMean_list<-lapply(EPMean_list,function(x) x[reord])
EPSD_list<-lapply(EPSD_list,function(x) x[reord])
CoVa_list<-lapply(CoVa_list,function(x) x[reord])


#We need to calculate the cell IDs for each slice, as the combined output only has a culmulative count of all cells in the volume
cellid<-function(rowNum){
	rowNum<-as.numeric(rowNum)
	slice<-outputRCT[rowNum,1]
	if (!slice == min(outputRCT[,1])){
		subtract<-max(which(outputRCT[,1]==slice-1))
		cellID<-outputRCT[rowNum,2]-outputRCT[subtract,2]
#the cell ID for the most ventral slice is just the cumulative cell ID
	}else{
		cellID<-outputRCT[rowNum,2]+1
	}
	returnList<-list("slice"=slice,"cellID"=cellID)
	return(returnList)
}

#Look at how the response (to the dots) change across the slice
outputS_list<-lapply(outputS_list,function(x) t(apply(x[,2:ncol(outputS_list[[1]])],1, function(x) x/mean(x))))
sliceMean<-t(sapply(outputS_list,colMeans))
sliceSD<-t(sapply(outputS_list,function(x) apply(x,2,sd)))

#Look at the ratio of the responses to dots and bars across slices (bar spans the whole visual field, whilst dot is only a small part)
ResponseRatio<-vector("list",length(outputS_list))
for (i in 1:length(outputS_list)) ResponseRatio[[i]]<-outputS_list[[i]][,1:(nEP/2)*2]/outputS_list[[i]][,1:(nEP/2)*2-1]
ResponseRatioMean<-t(sapply(ResponseRatio,colMeans))
ResponseRatioSD<-t(sapply(ResponseRatio,function(x) apply(x,2,sd)))



#Functions--------------------------------------------------------------------------------------------------------------
#1 	plot_cell(cell)			Plots the calcium trajectory of one cell
#2 	plot_dotselec()			Plots a histogram of the selectivity of dots vs bars for all sizes and directions
#3 	plot_global()			Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions 
#4 	plot_corr()				Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
#5 	plot_hist_resp()		Plots histogram of responses for all stimuli and a summary of the mean response and coefficient of variation
#6 	plot_sizeselec()		Plots a histogram of size selectivities
#7 	plot_cellloc(cell)		Plots the location of a cell in the tectum
#8 	plot_tuning(cell)		Plots a cell's response, dot vs bar selectivity, direction selectivity and location in the tectum
#9 	Plot_slice_mean()		Plots the mean response to the epochs in all slices
#10 Plot_ResponseRatio()	Plots the ratio of the dot to bar response across all slices


#1 Plots the calcium trajectory of one cell
plot_cell <- function(id,xmax=100000){
	param<-cellid(id)
	print(param)
	slice<-paste("_rslice",param$slice,sep="")
    par(mar=c(8,3,2,1))
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",param$cellID," ",folder,slice,sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix,id),xaxt='n',ann=F);
	axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2)
	rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(xmax,nEP*2),col=rgb(0,1,0,.05))
}

#2 Plots the histogram of the selectivity of dots vs bars for all sizes and directions
plot_dotselec<-function(){
	par(mfrow=c(2,5)); par(oma=c(0,0,4,0))

	for (i in 1:ncol(dotS_SF)){
		hist(dotS_SF[,i],breaks=seq(min(dotS_SF),max(dotS_SF),l=15),xlim=c(-1,1),ylim=c(0,600),xlab=paste(epochs[i*2],":",epochs[i*2-1],sep=""),main="",cex.lab=0.8)
	}
	mtext("dot to bar selectivity for each size and direction\n+1 = dot & -1 = bar",outer=T)
}

#Plot the histogram of DS for all sizes and bars/dots
plot_DS<-function(){
	    par(mfrow=c(2,5),oma=c(0,0,2,0));

	    for (i in 1:ncol(barDS_SF))
	        {
	        hist(barDS_SF[,i],breaks=seq(min(DSdf),max(DSdf),l=20),xlim=c(-1,1),ylim=c(0,280),xlab=paste(paste(EPlist[[2]][i]," - ",EPlist[[1]][i],sep="")),cex.lab=0.8,main="")
	        hist(dotDS_SF[,i],breaks=seq(min(DSdf),max(DSdf),l=20),xlim=c(-1,1),ylim=c(0,280),xlab=paste(paste(EPlist[[4]][i]," - ",EPlist[[3]][i],sep="")),cex.lab=0.8,main="")
	        }
	    mtext("Direction selectivity of dots and bars for each size\nvolumes", outer=T)
}


#3 Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions
plot_global<-function(P=19,S=1.5,thresh=0){
	par(mfrow=c(2,2)); par(oma=c(0,0,2,0));
	hist(dotDS,xlim=c(-1,1),ylim=c(0,200),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="Dot direction selectivity",xlab="90 degrees = +1 & 270 degrees = -1");
	hist(barDS,xlim=c(-1,1),ylim=c(0,200),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="Bar direction selectivity",xlab="90 degrees = +1 & 270 degrees = -1");
	hist(dotS90,xlim=c(-1,1),ylim=c(0,250),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (90 degrees)",xlab="dot = +1 & bar = -1");
	hist(dotS270,xlim=c(-1,1),ylim=c(0,250),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (270 degrees)",xlab="dot = +1 & bar = -1");
	mtext("volumes",outer=T)
	
#	for(i in 1:4){
#		globalsthresh[,i]<-globals[,i]>thresh|globals[,i]<(-thresh)
#		colPal<-colorRampPalette(colors=c("blue","white","red"))
#		cols<-colPal(50)[cut(globals[,i],breaks=50,labels=F)];
#		image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
#		points(x=cenfile[,2][globalsthresh[,i]],y=cenfile[,1][globalsthresh[,i]],col=cols,pch=P,cex=S)
#	}
}

#4 Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
plot_corr<-function(){
	par(mfrow=c(1,2));par(oma=c(0,0,2,0));
	plot(barDS ~ dotDS,xlim=c(-1,1),ylim=c(-1,1),main="bar DS vs dot DS")
	abline(DS.mod)
	plot(x=DS,y=dotS,xlim=c(-1,1),ylim=c(-1,1),main="direction selectivity vs dot selectivity")
	mtext("volumes",outer=T)
}

#5 Plots histogram of responses for all stimuli and a summary of the mean response and coefficient of variation
plot_hist_resp<-function(){

	par(mfrow=c(5,4));par(oma=c(0,0,2,0));
	for (i in 1:ncol(cellNorm)){
			hist(cellNorm[,i],xlab=epochs[i],xlim=c(0,0.1),ylim=c(0,750),breaks=seq(min(cellNorm),max(cellNorm),l=30),main="")
	}
	mtext("volumes",outer=T)

	X11();	lmat=matrix(c(1,2,3,4),ncol=2)
    layout(lmat)
    par(mar=c(6,3,3,3))
    par(oma=c(0,0,2,0))

	for (i in 1:length(EPsplit)){
		ymin<-EPMean_list[[i]]-EPSD_list[[i]]
		ymax<-EPMean_list[[i]]+EPSD_list[[i]]
		cent<-barplot(EPMean_list[[i]],xaxt="n",ylim=c(0,0.025))
		axis(side=1,labels=EPlist[[i]],cent,las=2,cex.axis=0.75)
	    arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,length=0.05,angle=90,code=3)
	}
	mtext("normalised mean response to epochs (volume)", outer=T)

	#plot coefficient of variation
	
	X11();  lmat=matrix(c(1,2,3,4),ncol=2)
    layout(lmat)
    par(mar=c(6,3,3,3))
    par(oma=c(0,0,2,0))
	
	for (i in 1:length(EPsplit)){
		cent<-barplot(CoVa_list[[i]],xaxt="n",ylim=c(0,0.65))
		axis(side=1,labels=EPlist[[i]],cent,las=2,cex.axis=0.75)
	}
	 mtext("normalised coefficient of variation (volume)", outer=T)
}

#6 Plots a histogram of size selectivities
plot_sizeselec<-function(){
	lmat=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T)
	layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5)); par(oma=c(0,0,2,0));
	minVal<-min(c(unlist(fivedeg),unlist(tendeg)))
	maxVal<-max(c(unlist(fivedeg),unlist(tendeg)))
	direction<-c("270","90")
	for (i in 1:3){
		for (j in 1:2) hist(fivedeg[[i]][,j],xlim=c(-1,1),ylim=c(0,400),breaks=seq(minVal,maxVal,l=20),xlab=paste(sizes[i+2]," degrees = -1 & 5 degrees = +1",sep=""),main=paste(direction[j],"degrees"))
	}
	mtext("5 degrees vs 17, 26 and 30 degrees (volume)",outer=T)

	 layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5)); par(oma=c(0,0,2,0));
	for (i in 1:3){
		for (j in 1:2) hist(tendeg[[i]][,j],xlim=c(-1,1),ylim=c(0,400),breaks=seq(minVal,maxVal,l=20),xlab=paste(sizes[i+2]," degrees = -1 & 10 degrees = +1",sep=""),main=paste(direction[j],"degrees"))
	}
	mtext("10 degrees vs 17, 26 and 30 degrees (volume)",outer=T)

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

		center<-barplot(as.numeric(Olist[[i]][cell,3:ncol(Olist[[i]])]),ylim=c(0,max(outputRCT[cell,3:ncol(outputRCT)])),main=gsub("_5",EPlist[[i]][1],replacement=""))
		axis(1,at=center,labels=EPlist[[i]],cex.axis=1,las=2)

	}	
		label1=c(paste(epochs[1:split*2],"/",epochs[1:split*2-1],sep=""),"DOT_90/BAR_90","DOT_270/BAR_270","DOT/BAR")
		label2=c(paste(EPlist[[2]],"-",EPlist[[1]],sep=""),paste(EPlist[[4]],"-",EPlist[[3]],sep=""),"dotDS","barDS","DS")
		col1=c(paste(rep("green4",ncol(dotS_SF)/2)),rep("blue",ncol(dotS_SF)/2),rep("red",2),rep("black",1))
		col2=c(paste(rep("green4",ncol(barDS_SF))),rep("blue",ncol(dotDS_SF)),rep("red",2),rep("black",1))

		dotchart(as.numeric(dotselec[cell,]),labels=label1,col=col1,pch=19,xlim=c(-1,1),main="dot/bar selectivity")
		dotchart(as.numeric(DSdf[cell,]),labels=label2,col=col2,pch=19,xlim=c(-1,1),main="Direction selectivity")

		plot_cell(cell)
		
		dc<-dev.cur()

		X11();
		plot_cellloc(cell);
		dev.set(dc)
}

plot_epoch_resp<-function(epoch,P=19,C=1.5){
	outputnorm<-apply(cellNorm,2,function(x) (x-min(x))/(max(x)-min(x)))
	image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
	points(x,y,col=rgb(1,0,0,outputnorm[,epoch]),pch=P,cex=C)
}

#9 Plots the mean response to the epochs in all slices
plot_slice_mean<-function(){
	ymin<-sliceMean-sliceSD
	ymax<-sliceMean+sliceSD
	
	par(mar=c(6,3,2,2))
	cent<-barplot(sliceMean,beside=T,ann=F,xaxt="n",main="mean response to epochs across imaging planes",ylim=c(0,max(ymax)))
	arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,code=3,length=0)
	axis(1,at=colMeans(cent),label=epochs,las=2,cex.axis=0.9)
}

#10 Plots the ratio of the dot to bar response across all slices
plot_response_ratio<-function(){
	ymin<-ResponseRatioMean-ResponseRatioSD
	ymax<-ResponseRatioMean+ResponseRatioSD

	par(mar=c(10,3,2,2))
	labs<-paste(epochs[1:(nEP/2)*2],"/",epochs[1:(nEP/2)*2-1],sep="")
	cent<-barplot(ResponseRatioMean,beside=T,ann=F,xaxt="n",main="ratio of dot to bar responses across imaging planes",ylim=c(0,max(ymax)))
	arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,code=3,length=0)
	axis(1,at=colMeans(cent),label=labs,las=2,cex.axis=0.8)
}

