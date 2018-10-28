martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");
#martincolorscale=c("green4","royalblue4","deeppink4","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

prefix_list=read.table("prefix_list_rep1")$V1

folder_list<-vector("list",length(prefix_list))
nifti_list<-vector("list",length(prefix_list))
outputRCT_list<-vector("list",length(prefix_list))
outputMC_list<-vector("list",length(prefix_list))
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
corvec_list<-vector("list",length(prefix_list))

for (i in 1:length(prefix_list)){
prefix<-prefix_list[i]
folder=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="");
epfile=paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="");

folder_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",prefix,sep="");
nifti_list[[i]]<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/rim_slice_1.nii",sep="");
outputRCT_list[[i]]<-read.table(paste(folder,prefix,"_output_rawCT.dat",sep=""));
noisethresh_list[[i]]<-read.table(paste(folder,prefix,"_noise_thresh.dat",sep=""));
outputMC_list[[i]]<-read.table(paste(folder,prefix,"_output_maxC.dat",sep=""))[noisethresh_list[[i]][,1],];
cenfile_list[[i]]<-read.table(paste(folder,prefix,"_cell_centers.dat",sep=""))[noisethresh_list[[i]][,1],];
avtab_list[[i]]<-read.table(paste(folder,prefix,"_mip.dat",sep=""));
dimvec_list[[i]]<-system(paste("nifti_tool -disp_hdr -infiles",nifti_list[[i]],"| grep \" dim \" | awk '{print $5, $6,$7}'"),intern=T);
dimvec_list[[i]]<-as.numeric(strsplit(dimvec_list[[i]]," ")[[1]]);
EPlength_list[[i]]<-read.table(paste(folder,prefix,"_EPtimeC.dat",sep=""))$V1
corvec_list[[i]]<-read.table(paste(folder,prefix,"_corvec.dat",sep=""))[noisethresh_list[[i]][,1],];
epochs<-read.table(paste(folder,"/",prefix,"_epochsC.dat",sep=""))$V1;

nx_list[[i]]=dimvec_list[[i]][1]
ny_list[[i]]=dimvec_list[[i]][2]
y_list[[i]]=cenfile_list[[i]]$V1+1
x_list[[i]]=cenfile_list[[i]]$V2+1

ep_list[[i]]<-read.table(epfile);

nEP<-nrow(ep_list[[i]])/2
ep_list[[i]]$V1<-floor((ep_list[[i]]$V1+90)*20);

#ord<-paste(gsub("\\D","",epochs),gsub("[^BD]","",epochs))
#ord<-order(gsub("3","03",ord))
ord<-order(epochs)
epochs<-epochs[ord]
outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])]<-outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])][,ord]

}

outputRCT<-do.call(rbind,outputRCT_list)
outputMC<-do.call(rbind,outputMC_list)
#output[[fish]][[epochs]]


#calculate the entropy of all the cells to all stimuli
outputRCT_list_ent<-lapply(outputRCT_list, function(output) t(apply(output[,2:ncol(output)],1,function(x) x/sum(x))))
outputRCT_list_ent<-lapply(outputRCT_list_ent, function(output) apply(output,1,function(x) sum(x*log2(1/x))))
max_ent<-log2(length(epochs))


rotate<-function(x) t(apply(x, 2, rev))

#extract global motion (SF15), dark and light local motion
index_GM<-grep("SF15",epochs)
index_LMD<-which("TRUE"==!grepl("LIGHTBAR",epochs)&grepl("BAR",epochs))
index_LML<-grep("LIGHTBAR",epochs)

epochs_GM<-epochs[index_GM]
epochs_LMD<-epochs[index_LMD]
epochs_LML<-epochs[index_LML]

outputRCT_GM<-outputRCT[,c(1,index_GM+1)]
outputRCT_LMD<-outputRCT[,c(1,index_LMD+1)]
outputRCT_LML<-outputRCT[,c(1,index_LML+1)]

outputMC_GM<-outputMC[,c(1,index_GM+1)]
outputMC_LMD<-outputMC[,c(1,index_LMD+1)]
outputMC_LML<-outputMC[,c(1,index_LML+1)]

#calculate selectivity indicies for global vs local motion
#NB we may need to normalise the global motion to eg the number of cycles of the grating....
dark_local_selec<-(outputRCT_LMD[,2:ncol(outputRCT_LMD)]-outputRCT_GM[,2:ncol(outputRCT_GM)])/(outputRCT_LMD[,2:ncol(outputRCT_LMD)]+outputRCT_GM[,2:ncol(outputRCT_GM)])
light_local_selec<-(outputRCT_LML[,2:ncol(outputRCT_LML)]-outputRCT_GM[,2:ncol(outputRCT_GM)])/(outputRCT_LML[,2:ncol(outputRCT_LML)]+outputRCT_GM[,2:ncol(outputRCT_GM)])
dark_max<-(outputMC_LMD[,2:ncol(outputMC_LMD)]-outputMC_GM[,2:ncol(outputMC_GM)])/(outputMC_LMD[,2:ncol(outputMC_LMD)]+outputMC_GM[,2:ncol(outputMC_GM)])
light_max<-(outputMC_LML[,2:ncol(outputMC_LML)]-outputMC_GM[,2:ncol(outputMC_GM)])/(outputMC_LML[,2:ncol(outputMC_LML)]+outputMC_GM[,2:ncol(outputMC_GM)])


#Functions--------------------------------------------------------------------------------------------------------------

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

plot_image<-function(cell,experiment,P=19,S=1.5){
	y=cenfile_list[[experiment]]$V1+1
	x=cenfile_list[[experiment]]$V2+1
	nx=dimvec_list[[experiment]][1]
	ny=dimvec_list[[experiment]][2]
	image(1:nx,1:ny,matrix(avtab_list[[experiment]]$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE);
	l<-outputRCT_list[[experiment]][,1]%in%(cell)
	points(x[l],y[l],col=martincolorscale[1],pch=P,cex=S)
			
}


#could maybe add in each tectum image...
plot_response_image<-function(EP,P=19,S=1.5){
	lmat<-matrix(1:length(prefix_list),ncol=length(prefix_list)/2)
	layout(lmat)

	for(i in 1:length(prefix_list)){
		nx<-nx_list[[i]]
		ny<-ny_list[[i]]
		avtab<-avtab_list[[i]]
		cenfile<-cenfile_list[[i]]
        colPal<-colorRampPalette(colors=c("blue","white","red"))
        cols<-colPal(50)[cut(outputRCT_list[[i]][,EP],breaks=50,labels=F)];
        image(1:nx,1:ny,matrix(avtab$V2,ncol=ny,nrow=nx),ylim=c(ny+0.5,0.5),col=gray.colors(100),xaxt='n',yaxt='n',ann=FALSE)
        points(x=cenfile[,2],y=cenfile[,1],col=cols,pch=P,cex=S)
     }
}




#9 plot the entropy of cells' response to all stimuli
plot_entropy<-function(){
	cols<-martincolorscale[1:length(prefix_list)]
	increment<-1/length(prefix_list)
	wex<-increment-(increment*0.1)
	loc<-(-0.45)
	boxplot(outputRCT_list_ent,las=2,cex.axis=0.8,boxfill=cols,ylim=c(0,max_ent_all),ylab="entropy")
	abline(h=max_ent_all,lty=3)
}

location_max<-function(){

	location_max3<-vector("list",length(prefix_list));
	location_max15<-vector("list",length(prefix_list));
	location_max30<-vector("list",length(prefix_list));
	epochs_W3 <-grep(("WIGGLY3_"), epochs, value = F, invert=F);
	epochs_W15 <-grep(("WIGGLY15_"), epochs, value = F, invert=F);
    epochs_W30 <-grep(("WIGGLY30_"), epochs, value = F, invert=F);

	for (i in 1:length(prefix_list)){
		location_max3[[i]]<-table(factor(epochs[epochs_W3[max.col(outputRCT_list[[i]][,epochs_W3+1])]],levels=epochs[epochs_W3]))
		location_max15[[i]]<-table(factor(epochs[epochs_W15[max.col(outputRCT_list[[i]][,epochs_W15+1])]],levels=epochs[epochs_W15]))
		location_max30[[i]]<-table(factor(epochs[epochs_W30[max.col(outputRCT_list[[i]][,epochs_W30+1])]],levels=epochs[epochs_W30]))
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

plot_local_selec<-function(){
	par(mfrow=c(2,2));
	minval<-min(c(min(dark_local_selec),min(light_local_selec)));
	maxval<-max(c(max(dark_local_selec),max(light_local_selec)));
	for(i in 1:ncol(dark_local_selec)) hist(dark_local_selec[,i],breaks=seq(minval,maxval,l=30),xlab=paste(epochs_LMD[i],"= +1 & ",epochs_GM[i],"= -1",sep=""),ylim=c(0,300),main="")
	X11(); par(mfrow=c(2,2));
	for(i in 1:ncol(light_local_selec)) hist(light_local_selec[,i],breaks=seq(minval,maxval,l=30),xlab=paste(epochs_LML[i],"= +1 & ",epochs_GM[i],"= -1",sep=""),ylim=c(0,300),main="")
	}

plot_local_selec_max<-function(){
	param_list<-vector("list",length=ncol(dark_max)+ncol(light_max))
	selec_list<-vector("list",length=2)
	par(mfrow=c(2,2),oma=c(4,0,4,0));
	minval<-min(c(min(dark_max),min(light_max)));
	maxval<-max(c(max(dark_max),max(light_max)));
	for(i in 1:ncol(dark_max)) param_list[[i]]<-hist(dark_max[,i],breaks=seq(minval,maxval,l=30),xlab="",ylim=c(0,300),main="",col=martincolorscale[i])
	mtext("global vs local motion selectivity (dark bar)",outer=T)
	par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
	plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
	legend("bottomleft",inset=c(0.01,0.01),legend=c("0ᵒ","180ᵒ","270ᵒ","90ᵒ"),xpd=TRUE,fill=martincolorscale[1:4],cex=1.3)
	legend("bottomright",inset=c(0.01,0.01),legend=c("1 = local selective","-1 = global selective"),xpd=TRUE,cex=1.3)
	
	X11(); par(mfrow=c(2,2),oma=c(4,0,4,0));
	for(i in 1:ncol(light_max)) param_list[[i+4]]<-hist(light_max[,i],breaks=seq(minval,maxval,l=30),xlab="",ylim=c(0,300),main="",col=martincolorscale[i])
	mtext("global vs local motion selectivity (light bar)",outer=T)
	par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
	plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
	legend("bottomleft",inset=c(0.01,0.01),legend=c("0ᵒ","180ᵒ","270ᵒ","90ᵒ"),xpd=TRUE,fill=martincolorscale[1:4],cex=1.3)
	legend("bottomright",inset=c(0.01,0.01),legend=c("1 = local selective","-1 = global selective"),xpd=TRUE,cex=1.3)
	
	selec_list[[1]]<-sapply(param_list, function(x) sum(x$counts[x$mid>=0.5]))
	selec_list[[2]]<-sapply(param_list, function(x) sum(x$counts[x$mid<=-0.5]))
	return(selec_list)
	}

low_ent_cells<-function(fish,amount){
	return(outputRCT_list[[fish]][outputRCT_list_ent[[fish]]%in%lapply(outputRCT_list_ent,sort)[[fish]][1:amount],1])
}

#from 0
cell_id<-function(X,experiment){
	return(outputRCT_list[[experiment]][X,1])
}

get_SINS<-function(experiment,xgt=180,ylt=200,ygt=70){
	y=cenfile_list[[experiment]]$V1+1
	x=cenfile_list[[experiment]]$V2+1
	id<-(cell_id(x>xgt & y<ylt & y>ygt,experiment))
	plot_image(id,experiment)
	return(id)
}

