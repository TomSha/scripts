require(plotrix)
require(mmand)

pulledresultsfolder=(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/pulled_results/",sep=""));
pulledoutput=read.table(paste(pulledresultsfolder,"/","all_centers.dat",sep=""));
pulledinc=read.table(paste(pulledresultsfolder,"/","pulled_inc.dat",sep=""))
pulleddencl=read.table(paste(pulledresultsfolder,"/","pulled_summary.dat",sep=""))
pulleddencl_sort=pulleddencl[order(pulleddencl$V1),];
#pulleddencl_inc=pulleddencl[pulledinc[,1],];
pulledoutput_inc=pulledoutput[pulledinc[,1],];
toolbox=paste("/home/meyer-lab/ToolBox_v3.6/")
prefix_list=read.table(paste(toolbox,"/scripts/prefix_list_rep1",sep=""))$V1

inc<-vector("list",length(prefix_list));
output<-vector("list",length(prefix_list));
dencl<-vector("list",length(prefix_list));
epochsC<-vector("list",length(prefix_list));
epochsnoC<-vector("list",length(prefix_list));
epochstime_list<-vector("list",length(prefix_list));
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

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i];
	resultsfolder<-(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep=""));
	timelogfolder<-(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep=""));
	bin_list[[i]]<-paste(resultsfolder,prefix,sep="");
	nifti_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");
	avtab_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_mip.dat",sep=""));

	epochsC[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_epochsC.dat",sep=""));
	ord<-order(epochsC[[i]]);
	epochsC[[i]]<-epochsC[[i]][,1][ord];
	epochsnoC[[i]]<-read.table(timelogfolder)$V2;
	nEP<-length(epochsnoC[[i]])/2;
	epochsnoC[[i]]<-epochsnoC[[i]][1:nEP*2];
	ord2<-order(epochsnoC[[i]]);
	epochstime_list[[i]]<-read.table(timelogfolder)$V1+90;
	epochstime_list[[i]]<-epochstime_list[[i]]*20; #multiple by Hz to get frame number

	output0[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_outputC.dat",sep=""));
	outputR0[[i]]<-read.table(paste(resultsfolder,prefix,"_output_rawC.dat",sep=""));
	V1<-outputR0[[i]][,1];
	outputR0[[i]]<-outputR0[[i]][,2:ncol(outputR0[[i]])][ord];
	outputR0[[i]]<-cbind(V1,outputR0[[i]]);
	outputRnoC[[i]]<-read.table(paste(resultsfolder,prefix,"_output_raw.dat",sep=""));
	outputRnoC[[i]]<-outputRnoC[[i]][,2:ncol(outputRnoC[[i]])][ord2];

	dencl[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_summary.dat",sep=""));
	dencl[[i]]<-dencl[[i]][order(dencl[[i]]$V1),];
	inc[[i]]<-read.table(paste(resultsfolder,prefix,"_inc.dat",sep=""));
	dencl[[i]]<-dencl[[i]][inc[[i]][,1],];
	cenfile_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_cell_centers.dat",sep=""));

	output[[i]]<-output0[[i]][output0[[i]]$V1 %in% dencl[[i]]$V1,];
	cells<-output[[i]][,1];
	output[[i]]<-output[[i]][,2:ncol(output[[i]])][ord];
	output[[i]]<-setNames(output[[i]],epochsC[[i]]);
	output[[i]]<-cbind(cells,output[[i]]);
	outputR[[i]]<-outputR0[[i]][outputR0[[i]]$V1 %in% dencl[[i]]$V1,];
}


##########################################################################

#create data frame summarising clustering data

maxn<-function(n) function(x) sort(x,decreasing = T)[n];
selec<-vector("list",length(prefix_list));

for (k in 1:length(prefix_list)){
	
	#calucluate max epoch value
	maxepochvalue<-apply(output[[k]][,2:ncol(output[[k]])],1,maxn(1));
	maxepochname<-colnames(output[[k]])[max.col(output[[k]][2:ncol(output[[k]])])+1]
	#calculate entropy
	normoutput[[k]]<-outputR[[k]][,2:ncol(outputR[[k]])]/rowSums(outputR[[k]][,2:ncol(outputR[[k]])])
	ent[[k]]<-apply(normoutput[[k]],1,function(x) sum(x*log(1/x)))


    no<-nrow(output[[k]])
	selec[[k]]<-data.frame(maxepochname=rep(NA,no),entropy=rep(NA,no),cellID=rep(NA,no),C1=rep(NA,no),C2=rep(NA,no))
	for (m in 1:nrow(output[[k]]))
	{
		prefix=prefix_list[k]
		selec[[k]]$maxepochname[m]	<- colnames(output[[k]])[which(maxepochvalue[m] == output[[k]][m,])];
		selec[[k]]$entropy[m]		<- ent[[k]][m];
		selec[[k]]$cellID[m] 		<- output[[k]][m,1];
		selec[[k]]$C1[m] 			<- dencl[[k]]$V6[m]+1

		if(any(pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m])){
			selec[[k]]$C2[m] <- pulleddencl_sort$V6[pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m]]+1
		} else {
			selec[[k]]$C2[m] <-"NULL"
		}
			
	}
}



#########################################################################
##FUNCTIONS##

#####print cells of a specific cluster

cells_in_cluster <-function(cluster){
	a <- vector("list",length(prefix_list))
	percent<-	"NULL"
	cell_num <- "NULL"
	for (j in 1:length(prefix_list)){
		a[[j]]<-selec[[j]][which(selec[[j]][,5]==cluster),]
		cell_num[j]<-nrow(a[[j]])						}
#	print(a)	
	print(cell_num)
	sum(as.integer(cell_num))
	}
#####


#####plots detrended single cell response

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
	ytop=rep(50000,nEP);

	system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID,bin));
    system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
	t<-read.table("out_dt.dat")$V1;
#	X11();
	plot(t,xaxt="n",type="l",ann=FALSE);
	axis(1,at=epochslab,labels=epochnames,las=2,cex.axis=0.7,tick=FALSE);
    axis(1,at=epochstime,labels=FALSE,tick=TRUE);
    rect(xleft=xleft,xright=xright,ybottom=ybot,ytop=ytop,col=rgb(0,0,1,.1));

}
#####


#####plot the detrended and normalised mean values for a cell

plot_dt_mean<-function(cell,exp){
	X11();
	par(mfrow=c(3,1));
	cell_ID<-which(output[[exp]]$cells==cell)
	plot_cell(cell+1,exp);
	barplot(as.matrix(outputR[[exp]][cell_ID,2:ncol(outputR[[exp]])]),xaxt="n",ann=FALSE);
	par(mar=c(8,4,0,2));
	barCenters<-barplot(as.matrix(output[[exp]][cell_ID,2:ncol(output[[exp]])]),xaxt="n",ann=FALSE);
	axis(1,at=barCenters,labels=epochsC[[exp]],las=2,cex.axis=0.7,line=2);
}
#####

###plot the mean raw response of all cells in a cluster across one fish
#plot_mean_cluster<-function(cluster, exp){
#	cell_ID  <- "NULL"
#	exp      <- selec[[experiment]][which(selec[[experiment]][,"C2"]==cluster),]
#	if (nrow(exp)>0){
#		for (q in 1:nrow(exp)){





#####plot the mean response of all cells in a cluster across all fish
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


#####plots the mean response of all cells a cluster for a given experiment (fish) 
#####cells have been normalised to their mean value
plot_dt_c<- function(cluster,experiment){

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
		
		cell_ID[q]<-exp[q,3]+1
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/read",cell_ID[q],bin));
		system(paste("/home/meyer-lab/ToolBox_v3.4/bin/detrend out_resp.dat out_dt.dat",dims[3],"20 5"));
		bin_cluster[[q]]<-read.table("out_dt.dat")$V1;
		bin_cluster[[q]]<-bin_cluster[[q]]/mean(bin_cluster[[q]])
		}
		
	a<-sapply(bin_cluster,unlist)
	mean_bin<-rowMeans(a)
#	X11()
#	plot(mean_bin,xaxt="n",type="l", col="red",ann=FALSE)
	plot(bin_cluster[[1]],xaxt="n",type="l", col="black",ann=F)
	title(main=paste(prefix,"    ",cluster),cex.main=0.75)
	axis(1,at=epochslab,labels=epochnames,las=2,cex.axis=0.55555,tick=FALSE);
	axis(1,at=epochstime,labels=FALSE,tick=TRUE)
	rect(xleft=xleft,xright=xright,ybottom=ybot,ytop=ytop,col=rgb(0,0,1,.1))
	colour=rgb(0,0,0,alpha=0.3)
	
	for (w in 2:nrow(exp)){
		lines(bin_cluster[[w]], col=colour)}
	lines(mean_bin,col="red")
}
}
#####


#####plot the mean response of all cells in all clusters for 3 fish 
plot_mean_cs<-function(exp1=0,exp2=0){

	include<-numeric(length=2)
	include<-c(exp1,exp2)
	print(include)
	for (cl in 1:max(pulleddencl$V6+1)){
		pdf(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/Plots/raw_traces/cluster",cl,"_exps",exp1,"-",exp2,".pdf",sep=""))
		par(mfrow=c(3,1))
		
		for (e in 1:length(include)){
			plot_dt_c(cl,include[e])

		}
		dev.off()
	}
}
#####


#####plots the average entropy of each cluster across all fish
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
	ymax<-max(Clmean,na.rm=T);
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


#####plots the entropy of all the cells in all of the fish
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


#####max cells per quadrants
location_max<-function(){

	location_max3<-vector("list",length(prefix_list));
	location_max15<-vector("list",length(prefix_list));
	location_max30<-vector("list",length(prefix_list));
	
	for (i in 1:length(prefix_list)){
		epochs_W3 <-grep(("WIGGLY3_"), epochsnoC[[i]], value = F, invert=F);
		epochs_W15 <-grep(("WIGGLY15_"), epochsnoC[[i]], value = F, invert=F);
        epochs_W30 <-grep(("WIGGLY30_"), epochsnoC[[i]], value = F, invert=F);
		outputRnoC[[i]]<-setNames(outputRnoC[[i]],epochsnoC[[i]]);

		location_max3[[i]]<-table(colnames(outputRnoC[[i]][epochs_W3])[max.col(outputRnoC[[i]][epochs_W3])]);
		location_max15[[i]]<-table(colnames(outputRnoC[[i]][epochs_W15])[max.col(outputRnoC[[i]][epochs_W15])]);
		location_max30[[i]]<-table(colnames(outputRnoC[[i]][epochs_W30])[max.col(outputRnoC[[i]][epochs_W30])]);

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

#####plot cmdscale for all cells according to second round of clustering

plot_cmd_all<-function(){

	outputcmd<-do.call(rbind,output);
	clustercmd<-do.call(rbind,selec);
	cl<-as.integer(clustercmd[,5])
	include<-!is.na(cl);
	outputcmd<-outputcmd[,2:ncol(outputcmd)];
	plot(cmdscale(dist(outputcmd[include,])),col=rainbow(max(cl[include]))[cl[include]],pch=19)
}
#####

#####plot cmdscale only for the centers from the first round of clustering
plot_cmd<-function(){
	cl<-pulleddencl_sort$V6[pulledinc[,1]]+1
	outputcmd<-pulledoutput_inc[,5:ncol(pulledoutput_inc)];
	plot(cmdscale(dist(outputcmd)),col=rainbow(max(cl))[cl],pch=19)
}
#####
plot_cmd_clusters<-function(){
	par(mfrow=c(3,3))
	outputcmd<-do.call(rbind,output);
	clustercmd<-do.call(rbind,selec);
	cl<-as.integer(clustercmd[,5])
	include<-!is.na(cl);
	outputcmd<-outputcmd[,2:ncol(outputcmd)];
	for (i in 1:max(na.omit(cl))){
	plot(cmdscale(dist(outputcmd[include,]))[cl[include]==i,],xlim=c(-4,4),ylim=c(-4,3),pch=19,ann=F);
	title(main=paste("cluster",i,sep=""));
	}

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
