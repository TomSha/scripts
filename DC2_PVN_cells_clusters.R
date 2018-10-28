require(plotrix)
require(mmand)
martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

pulledresultsfolder=(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/pulled_results/",sep=""));#CHANGE HERE
pulledoutput=read.table(paste(pulledresultsfolder,"/","all_centers.dat",sep=""));
pulledinc=read.table(paste(pulledresultsfolder,"/","pulled_inc.dat",sep=""))
pulleddencl=read.table(paste(pulledresultsfolder,"/","pulled_summary.dat",sep=""))
pulleddencl_sort=pulleddencl[order(pulleddencl$V1),];
#pulleddencl_inc=pulleddencl[pulledinc[,1],];
pulledoutput_inc=pulledoutput[pulledinc[,1],];
toolbox=paste("/home/meyer-lab/ToolBox_v3.6/")
prefix_list=read.table(paste(toolbox,"/scripts/prefix_list_rep1",sep=""))$V1 #CHANGE HERE

inc<-vector("list",length(prefix_list));
output<-vector("list",length(prefix_list));
dencl<-vector("list",length(prefix_list));
epochsC<-vector("list",length(prefix_list));
epochsnoC<-vector("list",length(prefix_list));
epochstime_list<-vector("list",length(prefix_list));
EPtimeC<-vector("list",length(prefix_list));
output0<-vector("list",length(prefix_list));
center_list<-vector("list",length(prefix_list));
bin_list<-vector("list",length(prefix_list));
outputRnoC<-vector("list",length(prefix_list));
outputR0<-vector("list",length(prefix_list));
outputR<-vector("list",length(prefix_list));
ent<-vector("list",length(prefix_list));
normoutput<-vector("list",length(prefix_list));
nifti_list<-vector("list",length(prefix_list));
cenfile_list<-vector("list",length(prefix_list));
avtab_list<-vector("list",length(prefix_list));
corvec_list<-vector("list",length(prefix_list));
noisethresh<-vector("list",length(prefix_list));
outnoremoval<-vector("list",length(prefix_list));

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i];
	resultsfolder<-(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep=""));
	timelogfolder<-(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep=""));
	bin_list[[i]]<-paste(resultsfolder,prefix,sep="");
	nifti_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");
	avtab_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_mip.dat",sep=""));
	noisethresh[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_noise_thresh.dat",sep=""))$V1	
	epochsC[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_epochsC.dat",sep=""));
	ord<-order(epochsC[[i]]);
	epochsC[[i]]<-epochsC[[i]][,1][ord];
	epochsnoC[[i]]<-read.table(timelogfolder)$V2;
	nEP<-length(epochsnoC[[i]])/2;
	epochsnoC[[i]]<-epochsnoC[[i]][1:nEP*2];
	ord2<-order(epochsnoC[[i]]);
	epochsnoC[[i]]<-epochsnoC[[i]][ord2]#new
	epochstime_list[[i]]<-read.table(timelogfolder)$V1+90;
	epochstime_list[[i]]<-epochstime_list[[i]]*20; #multiple by Hz to get frame number
	EPtimeC[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_EPtimeC.dat",sep=""))$V1
	EPtimeC[[i]]<-EPtimeC[[i]][ord]
	output0[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_outputC.dat",sep=""))[noisethresh[[i]],];
	outnoremoval[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_output_rawC.dat",sep=""))
	outputR0[[i]]<-read.table(paste(resultsfolder,prefix,"_output_rawCT.dat",sep=""))#[noisethresh[[i]],]
	outputR0[[i]][,2:ncol(outputR0[[i]])]<-outputR0[[i]][,2:ncol(outputR0[[i]])][,ord]
	outputRnoC[[i]]<-read.table(paste(resultsfolder,prefix,"_output_raw.dat",sep=""))[noisethresh[[i]],]#[,-2];
	outputRnoC[[i]][,2:ncol(outputRnoC[[i]])]<-outputRnoC[[i]][,2:ncol(outputRnoC[[i]])]#[ord2];

	dencl[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_summary.dat",sep=""));
	dencl[[i]]<-dencl[[i]][order(dencl[[i]]$V1),];
	inc[[i]]<-read.table(paste(resultsfolder,prefix,"_inc.dat",sep=""));
	dencl[[i]]<-dencl[[i]][inc[[i]][,1],];
	cenfile_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_cell_centers.dat",sep=""));

	output[[i]]<-output0[[i]][output0[[i]]$V1 %in% dencl[[i]]$V1,];
	output[[i]][,2:ncol(output[[i]])]<-output[[i]][,2:ncol(output[[i]])][ord];
	output[[i]]<-setNames(output[[i]],c("V1",as.matrix(epochsC[[i]])));
	outputR[[i]]<-outputR0[[i]][outputR0[[i]]$V1 %in% dencl[[i]]$V1,];
	outputR[[i]][,2:ncol(outputR[[i]])]<-outputR[[i]][,2:ncol(outputR[[i]])][ord];
	outputR[[i]]<-setNames(outputR[[i]],c("V1",as.matrix(epochsC[[i]])));
	outnoremoval[[i]][,2:ncol(outnoremoval[[i]])]<- outnoremoval[[i]][,2:ncol(outnoremoval[[i]])][ord]
	corvec_list[[i]]<-read.table(paste(resultsfolder,prefix,"_corvec.dat",sep=""))
}

#########################################################################
outputV<-vector("list",length(prefix_list));
outputVa<-vector("list",length(prefix_list));
outputVb<-vector("list",length(prefix_list));
emean<-vector("list",length(prefix_list));
emean2<-vector("list",length(prefix_list));

for (i in 1:length(prefix_list)){

	emean[[i]]<-colMeans(outputR[[i]][,2:ncol(outputR[[i]])]);
	emean2[[i]]<-colMeans(outputRnoC[[i]][,2:ncol(outputRnoC[[i]])]);
	outputV[[i]]<-outputR[[i]][,2:ncol(outputR[[i]])]/emean[[i]][col(outputR[[i]][,2:ncol(outputR[[i]])])];
	outputV[[i]]<-cbind(outputR[[i]][,1],outputV[[i]])

	outputVa[[i]]<-t(apply(outputR[[i]][,2:ncol(outputR[[i]])],1,function(x) (x/mean(x))))
	outputVa[[i]]<-outputVa[[i]]/emean[[i]][col(outputVa[[i]])];
	outputVa[[i]]<-cbind(outputR[[i]][,1],outputVa[[i]])

	outputVb[[i]]<-t(apply(outputR[[i]][,2:ncol(outputR[[i]])],1,function(x) (x-mean(x))))
	outputVb[[i]]<-outputVb[[i]]/emean[[i]][col(outputVb[[i]])];
	outputVb[[i]]<-cbind(outputR[[i]][,1],outputVb[[i]])

}

outputH<-lapply(outputR, function (x) t(apply(x[,2:ncol(outputR[[1]])],1,function(x) x/mean(x))))
for (i in 1:length(outputH)){
	for (j in 1:ncol(outputH[[i]])){
		outputH[[i]][,j]<-outputH[[i]][,j]/EPtimeC[[i]][j]
	}
}

alloutputH<-do.call(rbind,outputH)
#apply k means (k=2) over all the epochs and all fish
kfunc<-function(dat) kmeans(dat,centers=2)
km<-lapply(outputH,function(x) apply(x,2,kfunc))
cl1<-vector("list",length(prefix_list))
cl2<-vector("list",length(prefix_list))

for (i in 1:length(prefix_list)){
	for (j in 1:length(epochsC[[1]])){
		cl1[[i]][j]<-mean(outputH[[i]][,j][km[[i]][[j]]$cluster==1])
		cl2[[i]][j]<-mean(outputH[[i]][,j][km[[i]][[j]]$cluster==2])
	}
}
cl1<-do.call(rbind,cl1)
cl2<-do.call(rbind,cl2)


#look at outliers responses vs baseline responses
boxresults<-lapply(outputH,function(x) apply(x,2,function(x) boxplot(x,plot=F)))
outliermat<-matrix(rep(NA,length(prefix_list)*length(epochsC[[1]])),ncol=length(epochsC[[1]]),nrow=length(prefix_list))
medianmat<-matrix(rep(NA,length(prefix_list)*length(epochsC[[1]])),ncol=length(epochsC[[1]]),nrow=length(prefix_list))
for(i in 1:length(prefix_list)){
	for(j in 1:length(epochsC[[1]])){
		outliermat[i,j]<-mean(boxresults[[i]][[j]]$out[boxresults[[i]][[j]]$out>boxresults[[i]][[j]]$stats[5,1]])
		medianmat[i,j]<-boxresults[[i]][[j]]$stats[3,1]
	}
}

varvec<-t(sapply(outputH,function(x) apply(x,2,var)))



##########################################################################

#create data frame summarising clustering data

maxn<-function(n) function(x) sort(x,decreasing = T)[n];
selec<-vector("list",length(prefix_list));

for (k in 1:length(prefix_list)){
	
	#calucluate max epoch value
	maxepochvalue<-apply(outputR[[k]][,2:ncol(outputR[[k]])],1,maxn(1));
	maxepochname<-colnames(outputR[[k]])[max.col(outputR[[k]][2:ncol(outputR[[k]])])+1]
	#calculate entropy
	normoutput[[k]]<-outputR[[k]][,2:ncol(outputR[[k]])]/rowSums(outputR[[k]][,2:ncol(outputR[[k]])])
	ent[[k]]<-apply(normoutput[[k]],1,function(x) sum(x*log2(1/x)))


    no<-nrow(output[[k]])
	selec[[k]]<-data.frame(maxepochname=rep(NA,no),entropy=rep(NA,no),cellID=rep(NA,no),C1=rep(NA,no),C2=rep(NA,no),exp=rep(NA,no))
	for (m in 1:nrow(output[[k]]))
	{
		prefix=prefix_list[k]
		selec[[k]]$maxepochname[m]	<- colnames(outputR[[k]])[which(maxepochvalue[m] == outputR[[k]][m,])];
		selec[[k]]$entropy[m]		<- ent[[k]][m];
		selec[[k]]$cellID[m] 		<- output[[k]][m,1]+1
		selec[[k]]$C1[m] 			<- dencl[[k]]$V6[m]+1

		if(any(pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m])){
			selec[[k]]$C2[m] <- pulleddencl_sort$V6[pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m]]+1
		} else {
			selec[[k]]$C2[m] <-"NULL"
		}
		selec[[k]]$exp[m]			<-k
			
	}
}

cluster<-do.call(rbind,selec)
cluster<-as.integer(cluster[,5])

outputRall<-do.call(rbind,outputR);
outputRall<-outputRall[,2:ncol(outputRall)]
allcells<-do.call(rbind,selec)

selec3<-apply(outputRall,1,function(x) sum(x[grep("3_",epochsC[[1]])])/sum(x));				#FOR BARRAGE 2 (DOTS)
#selec3<-apply(outputRall,1,function(x) sum(x[grep("SF3",epochsC[[1]])])/sum(x));			#FOR BARRAGE 1
selec3<-cbind(selec3,allcells[,3],allcells[,6])
selec3<-selec3[order(selec3[,1],decreasing=T),]

selec15<-apply(outputRall,1,function(x) sum(x[grep("15_",epochsC[[1]])])/sum(x));
#selec15<-apply(outputRall,1,function(x) sum(x[grep("SF15",epochsC[[1]])])/sum(x));
selec15<-cbind(selec15,allcells[,3],allcells[,6])
selec15<-selec15[order(selec15[,1],decreasing=T),]

selec30<-apply(outputRall,1,function(x) sum(x[grep("30_",epochsC[[1]])])/sum(x));
#selec30<-apply(outputRall,1,function(x) sum(x[grep("SF30",epochsC[[1]])])/sum(x));
selec30<-cbind(selec30,allcells[,3],allcells[,6])
selec30<-selec30[order(selec30[,1],decreasing=T),]

DS_RL<-apply(outputRall,1,function(x) (sum(x[grep("_90",epochsC[[1]])])-sum(x[grep("_270",epochsC[[1]])]))/sum(x));
DS_RL<-cbind(DS_RL,allcells[,3],allcells[,6]);

DS_UD<-apply(outputRall,1,function(x) (sum(x[grep("_0",epochsC[[1]])])-sum(x[grep("_180",epochsC[[1]])]))/sum(x));
DS_UD<-cbind(DS_UD,allcells[,3],allcells[,6]);

halfnEP<-length(epochsC[[1]])/2
LDselec<-vector("list",length=halfnEP)
for (i in 1:halfnEP) LDselec[[i]]<-(outputRall[,i]-outputRall[,i+halfnEP])/(outputRall[,i]+outputRall[,i+halfnEP])

#DS<-DS[order(DS[,1],decreasing=T),]



#########################################################################
##FUNCTIONS##
#(0)	plot_means() - plots the mean response of all cells in all fish to each stimulus
#(1)	cells_in_cluster(cluster) - print cells of a specific cluster
#(2)	plot_cell(cell,exp) - plot detrended single cell response
#(3)	plot_cell_2(cell,exp) - plot detrended single cell response, raw mean, and normalised mean response.
#(4) 	plot_mean_cluster_all(cluster) - plot the mean raw response of all cells in a cluster across all fish
#(5) 	plot_cells(cluster, experiment) - plots the calcium trace of all cells in a cluster for a given experiment (fish) 
#(7) 	cluster_entropy() - plots the average entropy of each cluster across all fish
#(8) 	all_cells_entropy() - plots the entropy of all the cells in all of the fish
#(9) 	location_max() - plots the (normalised) number of cells for which that quadrant had the highest response
#(10)  	plot_cmd_all() - plot cmdscale for all cells according to second round of clustering
#(11)	plot_cmd() - plot cmdsclae on 2nd round clustering, cluster centers only
#(12) 	plot_cmd_clusters() - plots cmdscale for all cells with separate plots for each cluster

plot_interactions<-function(exp){
	plot(as.data.frame(outputH[[exp]]),xlim=c(0,0.04),ylim=c(0,0.04))
}

######(0) plot boxplot of epoch responses

plot_boxplot<-function(){
	cols<-martincolorscale[1:length(epochsC[[1]])]
	inc<-1/length(prefix_list)
	wex<-inc-(inc*0.1)
	loc<-(-0.45)
	par(mar=(c(10,3,1,2)));
	boxplot(x=outputH[[1]],las=2,cex.axis=0.8,boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1))

	for(i in 1:length(prefix_list)){
	res<-   		boxplot(x=outputH[[i]],add=T,boxwex=wex,at=1:length(epochsC[[1]])+loc,boxfill=cols,xaxt="n",yaxt="n")
		loc<-loc+inc
	}

}

plot_MedianAndOutlier<-function(){
	par(mar=c(8,3,2,2))
	cent<-barplot(outliermat,beside=T,main="median epoch response (colour) and mean of outlier responses")
	barplot(medianmat,beside=T,add=T,col=martincolorscale[1:7])
	axis(1,at=colMeans(cent),labels=epochsC[[1]],las=2,cex.axis=0.8)

	}


plot_kmean<-function(){
	
	sd1<-apply(cl1,2,sd)
	sd2<-apply(cl2,2,sd)
	M1<-colMeans(cl1)
	M2<-colMeans(cl2)

	ymin1<-sd1-M1
	ymax1<-sd1+M1
	ymin2<-sd2-M2
	ymax2<-sd2+M2
	
	ylim1<-min(c(ymin1,ymin2))
	ylim2<-max(c(ymax1,ymax2))

	xpos<-1:length(epochsC[[1]])

	par(mar=c(8,3,2,2))
	plot(M1,pch=19,xaxt="n",ann=F,ylim=c(ylim1,ylim2))
	arrows(x0=xpos,x1=xpos,y0=ymin1,y1=ymax1,length=0.03,angle=90,code=3)
	axis(1,labels=epochsC[[1]],at=1:length(epochsC[[1]]),las=2,cex.axis=0.8)
	points(colMeans(cl2),pch=19)
	arrows(x0=xpos,x1=xpos,y0=ymin2,y1=ymax2,length=0.03,angle=90,code=3)
}

#####(1) print cells of a specific cluster

cells_in_cluster <-function(cluster){
	a <- vector("list",length(prefix_list))
	percent<-	"NULL"
	cell_num <- "NULL"
	for (j in 1:length(prefix_list)){
		a[[j]]<-selec[[j]][which(selec[[j]][,5]==cluster),]
		cell_num[j]<-nrow(a[[j]])						}
	print(a)	
	print(cell_num)
	sum(as.integer(cell_num))
	}
#####



#####(2) plots detrended single cell response
plot_cell<-function(cell,exp){
	cell_ID <- cell
	prefix	<-prefix_list[exp]
	epochstime<-epochstime_list[[exp]]
	epochnames<-epochsnoC[[exp]]
	bin<-bin_list[[exp]]
	nifti<-nifti_list[[exp]]
    dims<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T);
    dims<-as.numeric(strsplit(dims," ")[[1]]);
    epochslab=colSums(matrix(epochstime,nrow=2))/2;
	nEP=length(epochstime)/2;
	xleft=epochstime[(1:nEP)*2-1];
	xright=epochstime[(1:nEP)*2];
	ybot=rep(0,nEP);
	ytop=rep(12000,nEP);

	system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID,bin));
    system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
	t<-read.table("out_dt.dat")$V1;
#	X11();
	par(mar=c(6,3,2,2))
	plot(t,xaxt="n",type="l",ann=FALSE);
	axis(1,at=epochslab,labels=epochnames,las=2,cex.axis=0.7,tick=FALSE,font=2);
    axis(1,at=epochstime,labels=FALSE,tick=TRUE);
    rect(xleft=xleft,xright=xright,ybottom=ybot,ytop=ytop,col=rgb(0,0,1,.1));

}
#####


#####(3) plot detrended single cell response, raw mean, and normalised mean response.
plot_cell_2<-function(cell,exp){
#	X11();
	par(mfrow=c(4,1), omi=c(0.5,0.3,0,0), plt=c(0.1,0.9,0,0.7));
#	par(mar=c(8,4,0,2));
	cell_ID<-which(output[[exp]]$V1==cell)-1
	plot_cell(cell,exp);
	barplot(as.matrix(outputR[[exp]][cell_ID,2:ncol(outputR[[exp]])]),xaxt="n",ann=FALSE);
	barplot(as.matrix(outputV[[exp]][cell_ID,2:ncol(outputV[[exp]])]),xaxt="n",ann=FALSE);
	barCenters<-barplot(as.matrix(output[[exp]][cell_ID,2:ncol(output[[exp]])]),xaxt="n",ann=FALSE);
	axis(1,at=barCenters,labels=epochsC[[exp]],las=2,cex.axis=0.7,line=2,font=2);
}
#####



#####(4) plot the mean raw response of all cells in a cluster across all fish
plot_mean_cluster_all<-function(cluster){
	center <- vector("list",length(prefix_list))
	for (i in 1:length(prefix_list)){
		center[[i]]<-outputR[[i]][which(selec[[i]][,5]==cluster),]
	}
	center<-do.call(rbind,center)
	normCenter<-t(apply(center[,2:ncol(center)],1,function(x) x/mean(x)))
	center_sd<-apply(normCenter,2,sd)
	y1=colMeans(normCenter)-center_sd
	y2=colMeans(normCenter)+center_sd
	cent<-barplot(colMeans(normCenter),xaxt="n")
	axis(side=1,labels=epochsC[[1]],at=cent,las=2,cex.axis=0.75)
	arrows(cent,y1,cent,y2,code=3,length=0)
	}
	
#####


#####(5) plots the calcium trace of all cells in a cluster for a given experiment (fish) 
#####cells have been normalised to their mean value
plot_cells<- function(cluster,experiment){

	cell_ID  <- "NULL"
	exp 	 <- selec[[experiment]][which(selec[[experiment]][,"C2"]==cluster),]	
		
	if (nrow(exp)>0){
	bin_cluster <- vector("list",nrow(exp))

	prefix<-prefix_list[experiment]
	epochstime<-epochstime_list[[experiment]]
	epochnames<-read.table(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep=""))$V2
	epochnames=epochnames[seq(1,length(epochnames),2)]
	bin<-bin_list[[experiment]]
	nifti<-nifti_list[[experiment]]
	dims<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dims<-as.numeric(strsplit(dims," ")[[1]])
	epochslab=colSums(matrix(epochstime,nrow=2))/2
	nEP=length(epochstime)/2;
	xleft=epochstime[(1:nEP)*2-1];
	xright=epochstime[(1:nEP)*2];
	ybot=rep(0,nEP);
	ytop=rep(50000,nEP);

	for (q in 1:nrow(exp)){
		
		cell_ID[q]<-exp[q,3]
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID[q],bin));
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
		bin_cluster[[q]]<-read.table("out_dt.dat")$V1;
		bin_cluster[[q]]<-bin_cluster[[q]]/mean(bin_cluster[[q]])
		}
		
	a<-sapply(bin_cluster,unlist)
	mean_bin<-rowMeans(a)
#	X11()
	par(mar=c(8,3,2,1))
	plot(mean_bin,xaxt="n",type="l", col="red",ann=FALSE)
#	plot(bin_cluster[[1]],xaxt="n",type="l", col="black",ann=F)
#	title(main=paste(prefix,"    ",cluster),cex.main=0.75)
	axis(1,at=epochslab,labels=epochnames,las=2,cex.axis=0.7,tick=FALSE,font=2);
	axis(1,at=epochstime,labels=FALSE,tick=TRUE)
	rect(xleft=xleft,xright=xright,ybottom=ybot,ytop=ytop,col=rgb(0,0,1,.1))
	colour=rgb(0,0,0,alpha=0.3)
	
#	for (w in 2:nrow(exp)){
#		lines(bin_cluster[[w]], col=colour)}
#	lines(mean_bin,col="red")
}
}
#####


#####(6) plot the mean response of all cells in all clusters for 3 fish 
plot_cells_2<-function(exp1=0,exp2=0){

	include<-numeric(length=2)
	include<-c(exp1,exp2)
	print(include)
	for (cl in 1:max(pulleddencl$V6+1)){
#		pdf(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/Plots/raw_traces/cluster",cl,"_exps",exp1,"-",exp2,".pdf",sep=""))
		par(mfrow=c(3,1))
		
		for (e in 1:length(include)){
			plot_cells(cl,include[e])

		}
#		dev.off()
	}
}
#####

plot_cluster_raw<-function(cluster,exp){
	center<-outputR[[exp]][which(selec[[exp]][,5]==cluster),2:ncol(outputR[[exp]])]
#	center<-t(apply(center,1,function(x) x/mean(x)))
	center_sd<-apply(center,2,sd)
	y1=colMeans(center)-center_sd
	y2=colMeans(center)+center_sd
#	par(omi=c(0.5,0.3,0,0))
	cent<-barplot(colMeans(center),xaxt="n")
#	axis(side=1,labels=epochsC[[1]],at=cent,las=2,cex.axis=0.75,font=2)
	arrows(cent,y1,cent,y2,code=3,length=0)
}

plot_cluster_vert<-function(cluster,exp){
	center<-outputV[[exp]][which(selec[[exp]][,5]==cluster),2:ncol(outputV[[exp]])]
#  	center<-t(apply(center,1,function(x) x/mean(x)))
	center_sd<-apply(center,2,sd)
	y1=colMeans(center)-center_sd
	y2=colMeans(center)+center_sd
#	par(omi=c(0.5,0.3,0,0))
	cent<-barplot(colMeans(center),xaxt="n")
#	axis(side=1,labels=epochsC[[1]],at=cent,las=2,cex.axis=0.75,font=2)
	arrows(cent,y1,cent,y2,code=3,length=0)
}

plot_cluster<-function(cluster,exp){
	center<-output[[exp]][which(selec[[exp]][,5]==cluster),2:ncol(output[[exp]])]
#	center<-t(apply(center,1,function(x) x/mean(x)))
	center_sd<-apply(center,2,sd)
	y1=colMeans(center)-center_sd
	y2=colMeans(center)+center_sd
	cent<-barplot(colMeans(center),xaxt="n")
#	par(omi=c(0.5,0.3,0,0))
	center<-do.call(rbind,center)
	axis(side=1,labels=epochsC[[1]],at=cent,las=2,cex.axis=0.75,font=2)
	arrows(cent,y1,cent,y2,code=3,length=0)
}

plot_cluster_all<-function(cluster,exp){
	par(mfrow=c(4,1),mar=c(6,3,2,1));
	plot_cells(cluster,exp);
	plot_cluster_raw(cluster,exp);
	plot_cluster_vert(cluster,exp);
	plot_cluster(cluster,exp);
}



#####(7) plots the average entropy of each cluster across all fish
cluster_entropy <-function(){
	exp <- vector("list",length(prefix_list));
	entropy <- vector("list",length(prefix_list));
	Clmax<-"NULL"

	for (i in 1:length(selec)){
		Clmax[i]<-max(as.numeric(selec[[i]][,5]),na.rm=T)
	}

	for (j in 1:length(prefix_list)){
		for (k in 1:max(Clmax)){
			exp[[j]]<-selec[[j]][which(selec[[j]][,5]==k),]
			entropy[[j]][k]<-mean(exp[[j]][,2])						
		}
		}
		
	entropy0<-t(sapply(entropy,unlist));
	Clmean<-colMeans(entropy0,na.rm=T);
	Clsd<-apply(entropy0,2,sd,na.rm=T);
	x1<-c(1:max(Clmax));
	y1<-Clmean-Clsd;
	y2<-Clmean+Clsd;
	ymax<-max(y2,na.rm=T);
	ymin<-min(y1,na.rm=T);
	X11();
	par(mar=c(6,4,4,2));
	plot(c(1:max(Clmax)),Clmean,xaxt="n",ann=F,ylim=c(ymin,ymax));
	arrows(x1,y1,x1,y2,length=0.05,angle=90,code=3);
	axis(1,at=c(1:max(Clmax)));
	mtext("cluster numer", side=1, line=3);
	mtext("entropy", side=2, line=3);
	}
#####


#####(8) plots the entropy of all the cells in all of the fish
all_cells_entropy<-function(){
	
	normoutputR0<-vector("list",length(prefix_list));
	allent<-vector("list",length(prefix_list));

	for (i in 1:length(outputR0)){
		normoutputR0[[i]]<-outputR0[[i]][,2:ncol(outputR0[[i]])]/rowSums(outputR0[[i]][,2:ncol(outputR0[[i]])]);
		allent[[i]]<-apply(normoutputR0[[i]],1,function(x) sum(x*log(1/x)));
		experiment<-sprintf(paste(prefix_list[i],"_","%03d",sep=""),outputR0[[i]][1:nrow(outputR0[[i]]),1]);
		maxEpoch<-as.character(epochsC[[i]][max.col(normoutputR0[[i]][,2:ncol(normoutputR0[[i]])])+1])
		allent[[i]]<-cbind(experiment,allent[[i]],maxEpoch);
		}

		unallent<-do.call(rbind,allent);
		unallent<-unallent[order(unallent[,2]),];
		X11();
		hist(as.numeric(unallent[,2]),freq=T,ann=F)
		mtext("entropy",side=1,line=3)
		mtext("frequency",side=2,line=3)

}
#####


#####(9) max cells per quadrants
location_max<-function(){

	location_max3<-vector("list",length(prefix_list));
	location_max15<-vector("list",length(prefix_list));
	location_max30<-vector("list",length(prefix_list));
	
	for (i in 1:length(prefix_list)){
		epochs_W3 <-grep(("WIGGLY3_"), epochsC[[i]], value = F, invert=F);
		epochs_W15 <-grep(("WIGGLY15_"), epochsC[[i]], value = F, invert=F);
        epochs_W30 <-grep(("WIGGLY30_"), epochsC[[i]], value = F, invert=F);
		out<-setNames(outputR0[[i]][,2:ncol(outputR0[[i]])],epochsC[[i]]);

		location_max3[[i]]<-table(factor(colnames(out[,epochs_W3])[max.col(out[,epochs_W3])],levels=epochsC[[i]][epochs_W3]));
		location_max15[[i]]<-table(factor(colnames(out[,epochs_W15])[max.col(out[,epochs_W15])],levels=epochsC[[i]][epochs_W15]));
		location_max30[[i]]<-table(factor(colnames(out[,epochs_W30])[max.col(out[,epochs_W30])],levels=epochsC[[i]][epochs_W30]));

		}

		location_max3<-t(sapply(location_max3,unlist));
		location_max15<-t(sapply(location_max15,unlist));
		location_max30<-t(sapply(location_max30,unlist));
		
		location_max3n<-t(apply(location_max3,1,function(x) x/mean(x)));
		location_max15n<-t(apply(location_max15,1,function(x) x/mean(x)));
		location_max30n<-t(apply(location_max30,1,function(x) x/mean(x)));
		
		max3_sd<-apply(location_max3n,2,sd);
		max15_sd<-apply(location_max15n,2,sd);
		max30_sd<-apply(location_max30n,2,sd);
		
		X11();
		par(mar=c(6,4,4,2))
		ymin<-min(colMeans(location_max3n)-max3_sd);
		ymax<-max(colMeans(location_max3n)+max3_sd);
		barCenters<-barplot(colMeans(location_max3n),axisnames=F, ylim=c(ymin,ymax));
		y1<-colMeans(location_max3n)-max3_sd;
		y2<-colMeans(location_max3n)+max3_sd;
		arrows(barCenters,y1,barCenters,y2,length=0.05,angle=90,code=3);
		mtext(text = colnames(location_max3n), side=1,at=barCenters, line=3);

		X11();
		ymin<-min(colMeans(location_max15n)-max15_sd);
		ymax<-max(colMeans(location_max15n)+max15_sd);
		barCenters<-barplot(colMeans(location_max15n), axisnames=F, ylim=c(ymin,ymax));
		y1<-colMeans(location_max15n)-max15_sd;
		y2<-colMeans(location_max15n)+max15_sd;
		arrows(barCenters,y1,barCenters,y2,length=0.05,angle=90,code=3);
		mtext(text = colnames(location_max15n), side=1,at=barCenters, line=3);

		X11();
		ymin<-min(colMeans(location_max30n)-max30_sd);
		ymax<-max(colMeans(location_max30n)+max30_sd);
		barCenters<-barplot(colMeans(location_max30n), axisnames=F,ylim=c(ymin,ymax));
		y1<-colMeans(location_max30n)-max30_sd;
		y2<-colMeans(location_max30n)+max30_sd;
		arrows(barCenters,y1,barCenters,y2,length=0.05,angle=90,code=3);
		mtext(text = colnames(location_max30n), side=1,at=barCenters, line=3);

}
#####

#####(10) plot cmdscale for all cells according to second round of clustering

plot_cmd_all<-function(){

	outputcmd<-do.call(rbind,output);
	clustercmd<-do.call(rbind,selec);
	cl<-as.integer(clustercmd[,5])
#remove clusters we don't include 
#	cl[cl==8|cl==14|cl==17|cl==18]<-NA
	include<-!is.na(cl);
	outputcmd<-outputcmd[,2:ncol(outputcmd)];
	plot(cmdscale(dist(outputcmd[include,])),col=rainbow(max(cl[include]))[cl[include]],pch=19)
}
#####

#####(11) plot cmdscale only for the centers from the first round of clustering
plot_cmd<-function(){
	cl<-pulleddencl_sort$V6[pulledinc[,1]]+1
#	cl[cl==8|cl==14|cl==17|cl==18]<-NA
	include<-!is.na(cl);
	outputcmd<-pulledoutput_inc[,5:ncol(pulledoutput_inc)];
	plot(cmdscale(dist(outputcmd)),col=rainbow(max(cl[include]))[cl[include]],pch=19)
}
#####(12) plots cmdscale for all cells with separate plots for each cluster
plot_cmd_clusters<-function(){
	par(mfrow=c(4,5))
	outputcmd<-do.call(rbind,output);
	clustercmd<-do.call(rbind,selec);
	cl<-as.integer(clustercmd[,5])
#	cl[cl==8|cl==14|cl==17|cl==18]<-NA
	include<-!is.na(cl);
	outputcmd<-outputcmd[,2:ncol(outputcmd)];
	for (i in 1:max(na.omit(cl))){
	plot(cmdscale(dist(outputcmd[include,]))[cl[include]==i,],xlim=c(-4,4),ylim=c(-4,4),pch=19,ann=F);
	title(main=paste("cluster",i,sep=""));
	}

}

### calculates within cluster variance
CV<-function(){
	alloutput<-do.call(rbind,output)
	alloutput<-alloutput[,2:ncol(alloutput)]
	cluster<-do.call(rbind,selec)
	cluster<-as.integer(cluster[,5])
	cluster[cluster==8|cluster==14|cluster==17|cluster==18]<-NA
	cluster[is.na(cluster)]<-0
	epochmeans<-matrix(nrow=max(cluster),ncol=ncol(alloutput))
	
	for (cl in 1:max(cluster)){
		epochmeans[cl,]<-colMeans(alloutput[cluster==cl,])
	}

	variance<-vector("list",max(cluster));
	for (i in 1:max(cluster)){
		if(sum(cluster==i)>0)
			variance[[i]]<-(alloutput[cluster==i,] - epochmeans[i])^2
		else
			variance[[i]]<-0
	}
	
	WCV<-colSums(do.call(rbind,variance))

	globalmean<-colMeans(alloutput[cluster,])
	
	BCV<-matrix(ncol=ncol(epochmeans),nrow=nrow(epochmeans))
	for (i in 1:max(cluster)){
			BCV[i,]<-((epochmeans[i,]-globalmean)^2)*nrow(do.call(rbind,variance))
	}
	BCV<-colSums(BCV,na.rm=T)
}

	
### plots stats about clusters and fish
stats<-function(){

#TAKE OUT THE CLUSTERS THAT WE DON'T INCLUDE!!!!!!!!!!
	clusterNum<-c(1:max(cluster,na.rm=T))	
	num_fish<-"NULL"
	mean_cells<-"NULL"
	sd_cells<-"NULL"
	variance<-"NULL"
   
   	cl<-lapply(selec,"[", ,5)

	for (i in 1:max(cluster,na.rm=T)){	
		num_fish[i]<-sum(lapply(cl,function(x) sum(x==i))==!0)
		mean_cells[i]<-round(as.numeric(mean(sapply(cl,function(x) sum(x==i))[!sapply(cl,function(x) sum(x==i))==0])))
		sd_cells[i]<-round(as.numeric(sd(sapply(cl,function(x) sum(x==i))[!sapply(cl,function(x) sum(x==i))==0])))

	}
	
	cluster_df<-data.frame(clusterNum,num_fish,mean_cells,sd_cells,variance)
	
	fish<-c(1:length(prefix_list))

	cells_clustered<-sapply(cl,function(x) sum(!x=="NULL"))#NEED TO TAKE OUT THE CLUSTERS WE DON'T INCLUDE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	percent<-round(sapply(cl,function(x) sum(!x=="NULL"))/sapply(output0,function(x) nrow(x)),2)
	clusters_in_fish<-round(sapply(cl,function(x) length(table(x)))-1)#/rep(max(cluster,na.rm=T),length(prefix_list)))

	fish_df<-data.frame(fish,cells_clustered,percent,clusters_in_fish)

}
################

plot_hist_resp<-function(X=c(1:ncol(alloutputH))){
	vec<-X
	X11();
	par(mfrow=c(5,4));
	for (i in 1:length(vec)){
		hist(alloutputH[,vec[i]],freq=F,main=epochsC[[1]][vec[i]],ylim=c(0,1),xlim=c(0,max(alloutputH)),breaks=seq(min(alloutputH),max(alloutputH),l=200),xlab="response")
	}

	meanvec<-lapply(outputH,function(x) colMeans(x))
	meanvec<-do.call(rbind,meanvec)
	sdvec<-lapply(outputH,function(x) apply(x,2,sd))
	sdvec<-do.call(rbind,sdvec)
	
	ymin<-meanvec-sdvec
	ymax<-meanvec+sdvec

	X11(); par(mar=c(8,3,2,1));
#	cent<-barplot(apply(alloutputH,2,mean),xaxt="n",main="mean response to each stimulus, normalised to epoch length")
	cent<-barplot(meanvec,beside=T,col=martincolorscale[1:length(prefix_list)],ann=F,xaxt="n",ylim=c(0,max(ymax)))
	arrows(x0=cent,x1=cent,y0=ymin,y1=ymax,angle=90,code=3,length=0.0)
	axis(side=1,labels=epochsC[[1]],at=colMeans(cent),las=2,cex.axis=0.75)
	
	CofV<-sdvec/meanvec
	X11(); par(mar=c(8,3,2,1));
#   cent<-barplot(apply(alloutputH,2,mean),xaxt="n",main="mean response to each stimulus, normalised to epoch length")
	cent<-barplot(CofV,beside=T,col=martincolorscale[1:length(prefix_list)],ann=F,xaxt="n",ylim=c(0,max(CofV)),main="coefficient of variation")
	axis(side=1,labels=epochsC[[1]],at=colMeans(cent),las=2,cex.axis=0.75)


	
}

plot_para<-function(){
	par(mfrow=c(3,1)); 
	hist(selec3[,1],freq=F,ylim=c(0,10),xlim=c(0,1),breaks=seq(0,1,l=100))
	hist(selec15[,1],freq=F,ylim=c(0,10),xlim=c(0,1),breaks=seq(0,1,l=100))
	hist(selec30[,1],freq=F,ylim=c(0,10),xlim=c(0,1),breaks=seq(0,1,l=100))
	X11();	hist(DS_UD[,1],breaks=seq(-1,1,l=100)) 
	X11();	plot(x=DS_UD[,1],y=allcells[,5])
	abline(v=-0.5,col="blue")
	abline(v=0.5,col="blue")

	X11();  hist(DS_RL[,1],breaks=seq(-1,1,l=100))
	X11();  plot(x=DS_RL[,1],y=allcells[,5],main="Direction selectivity per cluster: -1 = 270 degrees, +1 = 90 degrees",xlab="direction selectivity",ylab="cluster number")
	abline(v=-0.5,col="blue")
    abline(v=0.5,col="blue")

}
##################

plot_var<-function(){
	par(mar=c(7,5,1,1));
	barcenters<-barplot(varvec,xaxt="n",ylab="variance",beside=T);
	axis(1,at=colMeans(barcenters),labels=epochsC[[1]],las=2,cex.axis=0.7);
}
###################
rotate <- function(x) t(apply(x, 2, rev));

plot_heatmaps<-function(exp){
	X11();	image(t(as.matrix(outputR[[exp]][,2:ncol(outputR[[exp]])])),xaxt="n",yaxt="n",ann=F)#,col=colorRampPalette(colors=c("blue","white","red"))(100))
	X11(); 	image(t(as.matrix(output[[exp]][,2:ncol(output[[exp]])])),xaxt="n",yaxt="n",ann=F)#,col=colorRampPalette(colors=c("blue","white","red"))(100))
	
	cl<-as.numeric(selec[[exp]][,5][!is.na(as.numeric(selec[[exp]][,5]))])
	layout(matrix(1:2,ncol=2),widths=c(20,5))

	#image(y=1:length(cl),z=matrix(cl[order(cl)],nrow=1),col=rainbow(max(cl)),xaxt="n",yaxt="n",ann=F)

	axis(2,at=pos,labels=cl[order(cl)][pos],ylab=F)

	ordoutput<-as.matrix(output[[exp]][,2:ncol(output[[exp]])][!is.na(as.numeric(selec[[exp]][,5])),])
	ordoutput<-as.matrix(ordoutput[order(cl),])
	rangecol=range(ordoutput)
	pos<-floor(by(1:length(cl),cl[order(cl)],mean))
	
	image(x=1:ncol(ordoutput),y=1:length(cl),z=t(ordoutput[,2:ncol(ordoutput)]),xaxt='n',yaxt='n',xlab="",ylab="",frame=F)#,breaks=seq(rangecol[1],rangecol[2],length.out=101))#,col=colorRampPalette(colors=c("blue","white","red"))(100))
	image(y=1:length(cl),z,col=grey.colors(n=2,start=0,end=1),xaxt="n",yaxt="n",ann=F)

}	

#####

plot_corvec_hist<-function(exp){
	minbreak<-min(do.call(rbind,corvec_list))
	maxbreak<-max(do.call(rbind,corvec_list))
	hist(as.matrix(corvec_list[[exp]]),freq=F,xlim=c(0.45,0.75),ylim=c(0,23),col=rgb(0,0,1,0.3),breaks=seq(minbreak,maxbreak,l=20),main="",xlab="correlation")
}
######
plot_LDselec<-function(){
	par(mfrow=c(3,2));
	par(oma=c(0,0,4,0))
	for (i in 1:halfnEP){
		
		minVal<-min(unlist(LDselec))
		maxVal<-max(unlist(LDselec))
		hist(LDselec[[i]],xlim=c(-1,1),ylim=c(0,400),breaks=seq(minVal,maxVal,l=20),xlab=paste(paste(epochsC[[1]][i]," - ",epochsC[[1]][i+halfnEP],sep="")),main="")
	}
	mtext("selectivity for dark vs light dots\n+1 = dark dot & -1 = light dot",outer=T)
}






#####
Fscale <- function(x){
	        (x-min(x))/(max(x)-min(x))
}
#####

#####
HM<-function(exp,X,P=19,S=1){
	nifti<-nifti_list[[exp]]
	avtab<-avtab_list[[exp]]
	cenfile<-cenfile_list[[exp]]

	dimvec<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dimvec<-as.numeric(strsplit(dimvec," ")[[1]])
	nx<-dimvec[1]
	ny<-dimvec[2]
	x<-cenfile$V1+1
	y<-cenfile$V2+1


	    par(mfrow=c(2,2))
	    a<-grep(paste("WIGGLY",X,"_",sep=""),EPlist)+1
	    for(i in 1:length(a)){
		    image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=grey.colors(100),xaxt='n',yaxt='n',ann=FALSE);
		    title(main=EPlist[a[i]-1],cex.main=0.75)
		    points(x,y,col=rgb(1,0,0,Fscale(data0[,a[i]])),pch=P,cex=S);
	    }
}





#####
get_matrix <- function(experiment){
	prefix	<-prefix_list[experiment]
	epochstime<-epochstime_list[[experiment]]
	epochnames<-epochsnoC[[experiment]]
	bin<-bin_list[[experiment]]
	nifti<-nifti_list[[experiment]]
    dims<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T);
    dims<-as.numeric(strsplit(dims," ")[[1]]);
    epochslab=colSums(matrix(epochstime,nrow=2))/2;
	nEP=length(epochstime)/2;
	xleft=epochstime[(1:nEP)*2-1];
	xright=epochstime[(1:nEP)*2];
	mat<-matrix(NA,nrow=nrow(output0[[experiment]]),ncol=dims[[3]]);

	for(cell_ID in 1:nrow(output0[[experiment]])){
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID,bin));
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
		mat[cell_ID,]=read.table("out_dt.dat")$V1;
	}
	return(mat)
}
#####
gp_plot <- function(t,experiment){
	prefix	<-prefix_list[experiment]
	epochstime<-epochstime_list[[experiment]]
	epochnames<-epochsnoC[[experiment]]
    epochslab=colSums(matrix(epochstime,nrow=2))/2;
	nEP=length(epochstime)/2;
	xleft=epochstime[(1:nEP)*2-1];
	xright=epochstime[(1:nEP)*2];

	plot(t,xaxt="n",type="l",ann=FALSE);
	axis(1,at=epochslab,labels=epochnames,las=2,cex.axis=0.7,tick=FALSE);
    axis(1,at=epochstime,labels=FALSE,tick=TRUE);
    rect(xleft=xleft,xright=xright,ybottom=ybot,ytop=ytop,col=rgb(0,0,1,.1));
}

#####
plot_start_variance<-function(experiment){            
	       
	cell_ID  <- 1:nrow(outputRnoC[[experiment]]);

	if (length(cell_ID)>0){
	bin_exp <- vector("list",length(cell_ID));

	prefix<-prefix_list[experiment]
	epochstime<-epochstime_list[[experiment]]
	nEP<-length(epochstime)/2
	nifti<-nifti_list[[experiment]]
	dims<-system(paste("nifti_tool -disp_hdr -infiles",nifti,"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T)
	dims<-as.numeric(strsplit(dims," ")[[1]])

    for (q in 1:length(cell_ID)){
		bin<-bin_list[[experiment]];
        system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID[q],bin));	
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
        bin_exp[[q]]<-read.table("out_dt.dat")$V1;
			
    }
	}

#calculate variance before epochs begin and then variance in the epochs

	begin<-epochstime[seq(1,length(epochstime),2)];
	finish<-epochstime[seq(2,length(epochstime),2)];

	unbin_list<-t(sapply(bin_exp,unlist));
	unbin_list<-t(apply(unbin_list,1,function(x) x/mean(x)));
	
	unbin_list_start<-unbin_list[,1:begin[1]];

	for (i in 1:(length(begin)-1)){
	         unbin_list[,(finish[i]+1):(begin[i+1]-1)]<-NaN
			     }
	

	unbin_list_epochs<-unbin_list[,colSums(is.na(unbin_list)) != nrow(unbin_list)]
	unbin_list_epochs<-unbin_list_epochs[,-c(1:begin[1]-1)]

startVar<-apply(unbin_list_start,1,var);
epochVar<-apply(unbin_list_epochs,1,var);
thresh<-epochVar-startVar;
X11();
plot(thresh);
}
#####

