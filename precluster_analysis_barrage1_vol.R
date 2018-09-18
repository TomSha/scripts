martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

prefix="171025_PVN_barrage_F1_1_v10"

folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox/",prefix,sep="")
epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="")
niftifolder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/",sep="");

EPlength<-read.table(paste(folder,"_rslice10_EPtimeC.dat",sep=""))$V1

epfileC<-paste(folder,"_rslice10_epochsC.dat",sep="")
outputRCT<-read.table(paste(folder,"_combined_output_rawCT.dat",sep=""))
output<-read.table(paste(folder,"_combined_outputCT.dat",sep=""))
noisethresh<-read.table(paste(folder,"_noise_thresh.dat",sep=""))
cenfile<-read.table(paste(folder,"_cell_centers.dat",sep=""))[noisethresh[,1],]
#avtab<-read.table(paste(folder,"_mip.dat",sep=""))
#dimvec<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
#dimvec<-as.numeric(strsplit(dimvec," ")[[1]])
ep<-read.table(epfile)
ep$V1<-floor((ep$V1+90)*4);

epochsC<-read.table(epfileC)$V1;
ord<-order(epochsC)
epochsC<-epochsC[ord]

nEP<-length(epochsC)

outputRCT[,3:ncol(outputRCT)]<-outputRCT[,3:ncol(outputRCT)][,ord]
EPlength<-EPlength[ord]
for (i in 3:ncol(outputRCT)) outputRCT[,i]<-outputRCT[,i]/EPlength[i-2]

minSlice<-min(outputRCT[,1])
numSlices<-max(outputRCT[,1])-minSlice
nifti_list<-vector("list",length=numSlices);
outputS_list<-vector("list",length=numSlices);
noisethresh_list<-vector("list",length=numSlices);
cenfile_list<-vector("list",length=numSlices);
avtab_list<-vector("list",length=numSlices);
dimvec_list<-vector("list",length=numSlices);


for (i in 1:numSlices){
	slice<-i+(minSlice-1)
	nifti_list[[i]]<-paste(niftifolder,"/im_slice_",slice,"/rim_slice_",slice,".nii",sep="")
	outputS_list[[i]]<-read.table(paste(folder,"_rslice",slice,"_output_rawCT.dat",sep="")) 
	outputS_list[[i]]<-outputS_list[[i]][,2:ncol(outputS_list[[i]])][,ord]
	noisethresh_list[[i]]<-read.table(paste(folder,"_rslice",slice,"_noise_thresh.dat",sep=""))
	cenfile_list[[i]]<-read.table(paste(folder,"_rslice",slice,"_cell_centers.dat",sep=""))[noisethresh_list[[i]]$V1,]
	avtab_list[[i]]<-read.table(paste(folder,"_rslice",slice,"_mip.dat",sep=""))
	dimvec_list[[i]]<-system(paste("nifti_tool -disp_hdr -infiles",nifti_list[[i]],"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dimvec_list[[i]]<-as.numeric(strsplit(dimvec_list[[i]]," ")[[1]])
	}



#extract global motion (SF15), dark and light local motion
index_GM<-grep("SF15",epochsC)
index_LMD<-which("TRUE"==!grepl("LIGHTBAR",epochsC)&grepl("BAR",epochsC))
index_LML<-grep("LIGHTBAR",epochsC)

epochsC_GM<-epochsC[index_GM]
epochsC_LMD<-epochsC[index_LMD]
epochsC_LML<-epochsC[index_LML]

outputRCT_GM<-outputRCT[,c(1,2,index_GM+2)]
outputRCT_LMD<-outputRCT[,c(1,2,index_LMD+2)]
outputRCT_LML<-outputRCT[,c(1,2,index_LML+2)]

#calculate selectivity indicies for global vs local motion
#NB we may need to normalise the global motion to eg the number of cycles of the grating....
dark_local_selec<-(outputRCT_LMD[,3:ncol(outputRCT_LMD)]-outputRCT_GM[,3:ncol(outputRCT_GM)])/(outputRCT_LMD[,3:ncol(outputRCT_LMD)]+outputRCT_GM[,3:ncol(outputRCT_GM)])
light_local_selec<-(outputRCT_LML[,3:ncol(outputRCT_LML)]-outputRCT_GM[,3:ncol(outputRCT_GM)])/(outputRCT_LML[,3:ncol(outputRCT_LML)]+outputRCT_GM[,3:ncol(outputRCT_GM)])

#calculate the coefficient of variation 
EPMeans<-colMeans(outputRCT[,3:ncol(outputRCT)])
EPSD<-apply(outputRCT[,3:ncol(outputRCT)],2,sd)
EPCoV<-EPSD/EPMeans


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


#Functions-----------------------------
#1		plot_cell(cell)		Plots the calcium trajectory of one cell
#2		plot_local_selec()	plots a histogram of local vs global selectivity for global vs dark, and global vs light local stimuli
#3		plot_scatter()		plots a 2D scatter plot of local vs global motion
#4		plot_CoV()			plots the coefficient of variation for all stimuli
#5		plot_mean()			plots the mean and SD for all stimuli
#6 		plot_sd()			plots the sd for all stimuli

#1 Plots the calcium trajectory of one cell
plot_cell <- function(id,xmax=100000){
	param<-cellid(id)
	print(param)
	slice<-paste("_rslice",param$slice,sep="")
    par(mar=c(8,3,2,1))
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",param$cellID," ",paste(folder,slice,sep=""),sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix,id),xaxt='n',ann=F);
	axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2)
	rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(xmax,nEP*2),col=rgb(0,1,0,.05))
}

#2 plots a histogram of local vs global selectivity for global vs dark, and global vs light local stimuli
plot_local_selec<-function(){
	par(mfrow=c(2,2));
	minval<-min(c(min(dark_local_selec),min(light_local_selec)));
	maxval<-max(c(max(dark_local_selec),max(light_local_selec)));
	for(i in 1:ncol(dark_local_selec)) hist(dark_local_selec[,i],breaks=seq(minval,maxval,l=30),main=paste(epochsC_LMD[i],"-",epochsC_GM[i]))
	X11(); par(mfrow=c(2,2));
	for(i in 1:ncol(light_local_selec)) hist(light_local_selec[,i],breaks=seq(minval,maxval,l=30),main=paste(epochsC_LML[i],"-",epochsC_GM[i]))
	}


#3 plots scatter plot of global motion vs local motion
plot_scatter<-function(){
	maxval<-max(outputRCT[,3:ncol(outputRCT)])
	lmat<-matrix(seq(1,8),nrow=2)
	layout(lmat); par(mar=c(4,4,2,2))
	for (i in 3:ncol(outputRCT_GM)){
		plot(outputRCT_GM[,i],outputRCT_LMD[,i],xlab=epochsC_GM[i-2],ylab=epochsC_LMD[i-2],xlim=c(0,maxval),ylim=c(0,maxval))
		plot(outputRCT_GM[,i],outputRCT_LML[,i],xlab=epochsC_GM[i-2],ylab=epochsC_LML[i-2],xlim=c(0,maxval),ylim=c(0,maxval))
		}
}

#4 plots the coefficient of variation for all stimuli
plot_CoV<-function(){
	par(mar=c(8,3,2,2));
	cent<-barplot(EPCoV,xaxt="n");
	axis(1,at=cent,labels=epochsC,las=2)
}

#5 plots the mean and sd for all stimuli
plot_mean<-function(){
	par(mar=c(8,3,2,2));
	
	sdmax<-EPMeans+EPSD
	sdmin<-EPMeans-EPSD

	cent<-barplot(EPMeans,xaxt="n",ylim=c(0,max(sdmax)));
	arrows(x0=cent,y0=sdmin,x1=cent,y1=sdmax,code=3,length=0)
	axis(1,at=cent,labels=epochsC,las=2)
}

# plots sd for all stimuli
plot_sd<-function(){
	par(mar=c(8,3,2,2));
	 cent<-barplot(EPSD,xaxt="n")
	 axis(1,at=cent,labels=epochsC,las=2)
}

