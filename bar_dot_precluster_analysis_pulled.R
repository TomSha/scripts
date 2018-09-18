martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

prefix_list=read.table("prefix_list_bars_dots")$V1

folder_list<-vector("list",length(prefix_list))
nifti_list<-vector("list",length(prefix_list))
outputRCT_list<-vector("list",length(prefix_list))
noisethresh_list<-vector("list",length(prefix_list))
cenfile_list<-vector("list",length(prefix_list))
avtab_list<-vector("list",length(prefix_list))
dimvec_list<-vector("list",length(prefix_list))
nx_list<-vector("list",length(prefix_list))
ny_list<-vector("list",length(prefix_list))
y_list<-vector("list",length(prefix_list))
x_list<-vector("list",length(prefix_list))
ep_list<-vector("list",length(prefix_list))
EPlength_list<-vector("list",length(prefix_list))

for (i in 1:length(prefix_list)){
prefix<-prefix_list[i]
folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="");
epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="");

folder_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",prefix,sep="");
nifti_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");
outputRCT_list[[i]]<-read.table(paste(folder,prefix,"_output_rawCT.dat",sep=""));
noisethresh_list[[i]]<-read.table(paste(folder,prefix,"_noise_thresh.dat",sep=""));
cenfile_list[[i]]<-read.table(paste(folder,prefix,"_cell_centers.dat",sep=""))[noisethresh_list[[i]][,1],];
avtab_list[[i]]<-read.table(paste(folder,prefix,"_mip.dat",sep=""));
dimvec_list[[i]]<-system(paste("nifti_tool -disp_hdr -infiles",nifti_list[[i]],"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T);
dimvec_list[[i]]<-as.numeric(strsplit(dimvec_list[[i]]," ")[[1]]);
EPlength_list[[i]]<-read.table(paste(folder,prefix,"_EPtimeC.dat",sep=""))$V1

nx_list[[i]]=dimvec_list[[i]][1]
ny_list[[i]]=dimvec_list[[i]][2]
y_list[[i]]=cenfile_list[[i]]$V1+1
x_list[[i]]=cenfile_list[[i]]$V2+1

ep_list[[i]]<-read.table(epfile);

nEP<-nrow(ep_list[[i]])/2
ep_list[[i]]$V1<-floor((ep_list[[i]]$V1+90)*20);

epochs<-ep_list[[i]][,2][2:nEP*2]
ord<-order(paste(gsub("\\D","",epochs),gsub("[^BD]","",epochs)))
epochs<-epochs[ord]
outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])]<-outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])][,ord]
}

outputRCT<-do.call(rbind,outputRCT_list)


#minus 1 because we need to exclude the concentric
split<-(nEP-1)/2

#calculate selectivity for dots compared to bars
#dot-bar/total
dotS_SF<-(outputRCT[,2:ncol(outputRCT)][,1:split*2]-outputRCT[,2:ncol(outputRCT)][,1:split*2-1])/(outputRCT[,2:ncol(outputRCT)][,1:split*2]+outputRCT[,2:ncol(outputRCT)][,1:split*2-1])
dotS90<-rowSums(dotS_SF[,1:5])/5
dotS270<-rowSums(dotS_SF[,6:10])/5
dotS<-rowSums(dotS_SF[,1:10])/10
dotselec<-data.frame(dotS_SF,dotS90,dotS270,dotS)

#calculate direction selectivity (for bars and dots)
EPsplit<-c("BAR.270","BAR.90","DOT.270","DOT.90")
reord<-c(5,1,2,3,4)

EPlist<-vector("list",length=length(EPsplit))
Olist<-vector("list",length=length(EPsplit))
EPindex<-vector("list",length=length(EPsplit))
for (i in 1:length(EPsplit)){
	EPlist[[i]]<-epochs[grep(EPsplit[i],x=epochs)]
	Olist[[i]]<-outputRCT[,c(1,grep(EPsplit[i],x=epochs)+1)]
	EPindex[[i]]<-grep(pattern=EPsplit[i],epochs)
}
EPlist<-lapply(EPlist,function(x) x[reord])
Olist<-lapply(Olist,function(x) x[,c(1,(reord+1))])

#output[[fish]][[epochs]]
outputRCT_list_split<-lapply(outputRCT_list,function(output) lapply(EPindex, function(index) output[,c(index+1)]))

barDS_SF<-(Olist[[2]][,2:ncol(Olist[[2]])]-Olist[[1]][,2:ncol(Olist[[1]])])/(Olist[[2]][,2:ncol(Olist[[2]])]+Olist[[1]][,2:ncol(Olist[[1]])])
dotDS_SF<-(Olist[[4]][,2:ncol(Olist[[4]])]-Olist[[3]][,2:ncol(Olist[[3]])])/(Olist[[4]][,2:ncol(Olist[[4]])]+Olist[[3]][,2:ncol(Olist[[3]])])
barDS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("BAR_90",epochs)])-sum(x[grep("BAR_270",epochs)]))/sum(x[grep("BAR",epochs)]))
dotDS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("DOT_90",epochs)])-sum(x[grep("DOT_270",epochs)]))/sum(x[grep("DOT",epochs)]))
DS<-apply(outputRCT[,2:ncol(outputRCT)],1,function(x) (sum(x[grep("90",epochs)])-sum(x[grep("270",epochs)]))/sum(x))
DSdf<-data.frame(barDS_SF,dotDS_SF,barDS,dotDS,DS)

DS.mod<-lm(barDS ~ dotDS)

globals<-data.frame(dotDS,barDS,dotS90,dotS270);

#look at selectivity comparing different size dots. Compare 5 vs 17,26 and 30 degrees. Compare 10 vs 17,26,30 degrees
sizes<-as.integer(substr(EPlist[[3]],9,10))
sizeList<-vector("list",length=length(sizes))
fivedeg<-vector("list",length=3)
tendeg<-vector("list",length=3)
for (i in 1:length(sizeList)) sizeList[[i]]<-outputRCT[,2:ncol(outputRCT)][grepl("DOT",epochs)&grepl(sizes[i],epochs)]
for (i in 3:length(sizeList)) fivedeg[[i-2]]<-(sizeList[[1]]-sizeList[[i]])/(sizeList[[1]]+sizeList[[i]])
for (i in 3:length(sizeList)) tendeg[[i-2]]<-(sizeList[[2]]-sizeList[[i]])/(sizeList[[2]]+sizeList[[i]])


#look at the mean responses across epochs
EPlength<-do.call(rbind,EPlength_list)
#normalise across cells to mean and then across epochs to epoch length
cellNorm_list<-sapply(outputRCT_list, function(x) t(apply(x[,2:ncol(outputRCT_list[[1]])],1,function(x) x/mean(x))))
for (i in 1:length(cellNorm_list)){
	for (j in 1:ncol(cellNorm_list[[i]])){
		cellNorm_list[[i]][,j]<-cellNorm_list[[i]][,j]/EPlength_list[[i]][j]
	}
}

EPMeans<-t(sapply(as.matrix(cellNorm_list),colMeans))
#applies sd over epochs for each fish
EPSD<-t(sapply(cellNorm_list, function(x) apply(x,2,sd)))

EPMean_list<-vector("list",length=length(EPsplit))
EPSD_list<-vector("list",length=length(EPsplit))
CoVa_list<-vector("list",length=length(EPsplit))
for (i in 1:length(EPsplit)){
	EPMean_list[[i]]<-EPMeans[,grep(EPsplit[i],x=epochs)]
	EPSD_list[[i]]<-EPSD[,grep(EPsplit[i],x=epochs)]
	CoVa_list[[i]]<-EPSD_list[[i]]/EPMean_list[[i]]
}

EPMean_list<-lapply(EPMean_list,function(x) x[,reord])
EPSD_list<-lapply(EPSD_list,function(x) x[,reord])
CoVa_list<-lapply(CoVa_list,function(x) x[,reord])

#Calculate the entropy of all the cells response to different sizes 
outputRCT_list_ent<-lapply(outputRCT_list_split, function(output) lapply(output, function(output) output/rowSums(output)))
outputRCT_list_ent<-lapply(outputRCT_list_ent, function(output) lapply(output,function(output) apply(output,1,function(x) sum(x*log2(1/x)))
max_ent<-log2(length(EPlist[[1]]))


#Functions--------------------------------------------------------------------------------------------------------------
#1 plot_cell(cell,exp)  Plots the calcium trajectory of one cell
#2 plot_dotselec()      Plots a histogram of the selectivity of dots vs bars for all sizes and directions
#3 plot_DS()			Plots the histogram of DS for all sizes and bars/dots
#4 plot_global()        Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions 
#5 plot_corr()          Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
#6 plot_hist_resp()     Plots histogram of responses for all stimuli and a summary of the mean response and coefficient of variation
#7 plot_sizeselec()     Plots a histogram of size selectivities

#1 Plots the calcium trajectory of one cell
plot_cell <- function(id,exp,xmax=100000){
	ep<-ep_list[[exp]]
    par(mar=c(8,3,2,1))
    system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder_list[[exp]],sep=""));
    t<-read.table("out_resp.dat")$V1;
    system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
    t_dt<-read.table("out_dt.dat")$V1;
    plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix_list[exp],id),xaxt='n',ann=F);
    axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2)
    rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(xmax,nEP*2),col=rgb(0,1,0,.05))
}

#2 Plots the histogram of the selectivity of dots vs bars for all sizes and directions
plot_dotselec<-function(){
    par(mfrow=c(2,5));par(oma=c(0,0,4,0));

	i=5; hist(dotS_SF[,i],breaks=seq(min(dotselec),max(dotselec),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(epochs[i*2],":",epochs[i*2-1],sep=""),main="",cex.lab=0.8)
    for (i in 1:4){
        hist(dotS_SF[,i],breaks=seq(min(dotselec),max(dotselec),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(epochs[i*2],":",epochs[i*2-1],sep=""),main="",cex.lab=0.8)
	    }
	i=10; hist(dotS_SF[,i],breaks=seq(min(dotselec),max(dotselec),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(epochs[i*2],":",epochs[i*2-1],sep=""),main="",cex.lab=0.8)
	for (i in 6:9){
		        hist(dotS_SF[,i],breaks=seq(min(dotselec),max(dotselec),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(epochs[i*2],":",epochs[i*2-1],sep=""),main="",cex.lab=0.8)
				  }
	mtext("dot to bar selectivity for each size and direction\n+1 = dot & -1 = bar", outer=T)
}

#3 Plots the histogram of DS for all sizes and bars/dots
plot_DS<-function(){
	par(mfrow=c(2,5),oma=c(0,0,2,0));

	for (i in 1:ncol(barDS_SF))
		{
	    hist(barDS_SF[,i],breaks=seq(min(DSdf),max(DSdf),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(paste(EPlist[[2]][i]," - ",EPlist[[1]][i],sep="")),cex.lab=0.8,main="")
		hist(dotDS_SF[,i],breaks=seq(min(DSdf),max(DSdf),l=20),xlim=c(-1,1),ylim=c(0,220),xlab=paste(paste(EPlist[[4]][i]," - ",EPlist[[3]][i],sep="")),cex.lab=0.8,main="")
	    }
	mtext("Direction selectivity of dots and bars for each size", outer=T)
}


#4 Plots a histogram of direction selectivity of bars and dots and dot selectivity for both directions
plot_global<-function(P=19,S=1.5,thresh=0){
    par(mfrow=c(2,2));
    hist(dotDS,xlim=c(-1,1),ylim=c(0,150),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="Dot direction selectivity",xlab="90 degrees = +1 & 270 degrees = -1");
    hist(barDS,xlim=c(-1,1),ylim=c(0,150),breaks=seq(min(c(dotDS,barDS)),max(c(dotDS,barDS)),l=30),main="Bar direction selectivity",xlab="90 degrees = +1 & 270 degrees = -1");
    hist(dotS90,xlim=c(-1,1),ylim=c(0,150),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (90 degrees)",xlab="dot = +1 & bar = -1");
    hist(dotS270,xlim=c(-1,1),ylim=c(0,150),breaks=seq(min(c(dotS90,dotS270)),max(c(dotS90,dotS270)),l=30),main="dot selectivity (270 degrees)",xlab="dot = +1 & bar = -1");

#could maybe add in each tectum image...
#	for(i in 1:4){
#        globalsthresh[,i]<-globals[,i]>thresh|globals[,i]<(-thresh)
#        colPal<-colorRampPalette(colors=c("blue","white","red"))
#        cols<-colPal(50)[cut(globals[,i],breaks=50,labels=F)];
#        image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
#        points(x=cenfile[,2][globalsthresh[,i]],y=cenfile[,1][globalsthresh[,i]],col=cols,pch=P,cex=S)
#	     }
}

#5 Plots dot vs bar direction selectivity and direction selectivity vs dot selectivity
plot_corr<-function(){
	par(mfrow=c(1,2));
    plot(barDS ~ dotDS,xlim=c(-1,1),ylim=c(-1,1),main="bar DS vs dot DS")
    abline(DS.mod)
    plot(x=DS,y=dotS,xlim=c(-1,1),ylim=c(-1,1),main="direction selectivity vs dot selectivity")
}

#6 Plots histogram of responses for epochs, all cells are normalised to their mean response, then each epoch normalised to the length of the epoch.  
plot_hist_resp<-function(){
	
	cellNorm<-do.call(rbind,cellNorm_list)		
	par(mfrow=c(5,4));
	for (i in 1:ncol(cellNorm)){
		hist(cellNorm[,i],main=epochs[i],xlim=c(0,0.1),ylim=c(0,500),breaks=seq(min(cellNorm[,1:ncol(cellNorm)]),max(cellNorm[,1:ncol(cellNorm)]),l=70))
}
	#plot mean of means 
	X11()
	lmat=matrix(c(1,2,3,4),ncol=2)
	layout(lmat)
	par(mar=c(6,3,3,3))
	par(oma=c(0,0,2,0))
	for (i in 1:length(EPsplit)) {
		EPM<-colMeans(EPMean_list[[i]])
		EPsd<-apply(EPMean_list[[i]],2,sd)
		ymin<-EPM-EPsd
		ymax<-EPM+EPsd
		cent<-barplot(EPM,xaxt="n",ylim=c(0,0.016))
		axis(side=1,labels=EPlist[[i]],cent,las=2,cex.axis=0.75)
		arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,length=0.05,angle=90,code=3)
	}
	mtext("normalised mean response to epochs",outer=T)

	#plot means of coefficient of variation
	X11()
	lmat=matrix(c(1,2,3,4),ncol=2)
	layout(lmat)
	par(mar=c(6,3,3,3))
	par(oma=c(0,0,2,0))
	for (i in 1:length(EPsplit)){
		CoVa<-colMeans(CoVa_list[[i]])
		CoVaSD<-apply(CoVa_list[[i]],2,sd)
		ymin<-CoVa-CoVaSD
		ymax<-CoVa+CoVaSD
		cent<-barplot(CoVa,xaxt="n",ylim=c(0,1.2))
		axis(side=1,labels=EPlist[[i]],cent,las=2,cex.axis=0.75)
		arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,length=0.05,angle=90,code=3)
	}
	mtext("mean of cofficient of variation",outer=T)
}

#7 Plots a histogram of size selectivities
plot_sizeselec<-function(){
	lmat=matrix(c(1,2,3,4,5,6),ncol=2,byrow=T)
	layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5));	par(oma=c(0,0,2,0)) 
	minVal<-min(c(unlist(fivedeg),unlist(tendeg)))
	maxVal<-max(c(unlist(fivedeg),unlist(tendeg)))
	direction<-c("270","90")
	for (i in 1:3){
		for (j in 1:2) hist(fivedeg[[i]][,j],xlim=c(-1,1),ylim=c(0,270),breaks=seq(minVal,maxVal,l=20),xlab=paste(sizes[i+2]," degrees = -1 & 5 degrees = +1 ",sep=""),main=paste(direction[j],"degrees"))
			    }
	mtext("5 degrees vs 17, 26 and 30 degrees",outer=T)
#	X11();
	layout(lmat)
	par(mar=c(6.5,2.5,2.5,2.5));	par(oma=c(0,0,2,0))
	for (i in 1:3){
		for (j in 1:2) hist(tendeg[[i]][,j],xlim=c(-1,1),ylim=c(0,270),breaks=seq(minVal,maxVal,l=20),xlab=paste(sizes[i+2]," degrees = -1 & 10 degrees = +1 ",sep=""),main=paste(direction[j],"degrees"))
			    }
	mtext("10 degrees vs 17, 26 and 30 degrees",outer=T)
}

#8 plot the entropy of cells' response across different sizes
plot_entropy<-function(){
	cols<-martincolorscale[1:length(EPsplit)]
	increment<-1/length(prefix_list)
	wex<-increment-(increment*0.1)
	loc<-(-0.45)
	boxplot(x=outputRCT_list_ent[[1]],las=2,cex.axis=0.8,boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),ylim=c(0,max_ent),ylab="entropy")
	abline(h=max_ent,lty=3)
	
	for(i in 1:length(prefix_list)){
		boxplot(lapply(outputRCT_list_ent, `[[`, i),add=T,boxwex=wex,at=1:length(EPsplit)+loc,xaxt="n",yaxt="n",boxfill=cols)
		loc<-loc+increment
	}
	legend("bottomleft",inset=0.05,legend=EPsplit,fill=martincolorscale[1:length(EPsplit)])

}



