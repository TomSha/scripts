martincolorscale=c("#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");
#martincolorscale=c("green4","royalblue4","deeppink4","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

prefix_list=read.table("prefix_list_dots")$V1

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
epochmod<-gsub("3_","03_",epochs)
ord<-order(paste(gsub("\\D","",epochmod),gsub("[^BD]","",epochmod)))
outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])]<-outputRCT_list[[i]][,2:ncol(outputRCT_list[[i]])][,ord]
}
epochmod<-epochmod[ord]
epochs<-epochs[ord]
outputRCT<-do.call(rbind,outputRCT_list)
split<-(nEP-1)/2

#calculate the entropy of all the cells to all stimuli
outputRCT_list_ent_all<-lapply(outputRCT_list, function(output) t(apply(output[,2:ncol(output)],1,function(x) x/sum(x))))
outputRCT_list_ent_all<-lapply(outputRCT_list_ent_all, function(output) apply(output,1,function(x) sum(x*log2(1/x))))
max_ent_all<-log2(length(epochs))

darkselec<-(outputRCT[,2:ncol(outputRCT)][,(1:split)*2]-outputRCT[,2:ncol(outputRCT)][,(1:split)*2-1])/(outputRCT[,2:ncol(outputRCT)][,(1:split)*2]+outputRCT[,2:ncol(outputRCT)][,(1:split)*2-1])
darkselec_mean<-apply(darkselec,1,mean)

#Functions--------------------------------------------------------------------------------------------------------------
#1 plot_cell(cell,exp)  Plots the calcium trajectory of one cell
#2 plot the entropy of cells' response to all stimuli

#1 Plots the calcium trajectory of one cell
plot_cell <- function(id,exp,xmax=100000){
	ep<-ep_list[[exp]]
    par(mar=c(8,3,2,1))
    system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder_list[[exp]],sep=""));
    t<-read.table("out_resp.dat")$V1;
    system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
    t_dt<-read.table("out_dt.dat")$V1;
    plot((t_dt),type='l',ylim=(c(0,max(t_dt))),main=paste(prefix_list[exp],id),xaxt='n',ann=F);
    axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2,cex.axis=0.8)
    rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(xmax,nEP*2),col=rgb(0,1,0,.05))
}

#2 plot the entropy of cells' response to all stimuli
plot_entropy_all<-function(){
	cols<-martincolorscale[1:length(prefix_list)]
	increment<-1/length(prefix_list)
	wex<-increment-(increment*0.1)
	loc<-(-0.45)
	boxplot(outputRCT_list_ent_all,las=2,cex.axis=0.8,boxfill=cols,ylim=c(0,max_ent_all),ylab="entropy")
	abline(h=max_ent_all,lty=3)
}

plot_darkselec<-function(){
	param_list<-vector("list",length=ncol(darkselec))
	selec_list<-vector("list",length=2)

#	dev.new(width=20,height=10)
	lmat<-matrix(1:6,ncol=3,byrow=F)
	layout(lmat)
    par(oma=c(6,0,4,0));
	param_list<-vector("list",length=ncol(darkselec))
	selec_list<-vector("list",length=2)
	for(i in 1:6){
		if(i%%2){
			param_list[[i]]<-hist(darkselec[,i],breaks=seq(min(darkselec),max(darkselec),l=20),xlim=c(-1,1),ylim=c(0,400),xlab="",main=paste(paste(substr(epochmod[i*2-1],start=12,stop=13),"ᵒ",sep="")),cex.main=1.5,,cex.lab=1,col=martincolorscale[1])
		}
		else
			param_list[[i]]<-hist(darkselec[,i],breaks=seq(min(darkselec),max(darkselec),l=20),xlim=c(-1,1),ylim=c(0,400),xlab="",cex.main=1.5,,cex.lab=1,col=martincolorscale[2],main="")
	}
	mtext("light vs dark dot selectivity", outer=T,cex=2)

	par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
	plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
	legend("bottomleft",inset=c(0.01,0.01),legend=c("270ᵒ direction","90ᵒ direction"),xpd=TRUE,fill=martincolorscale[1:2],cex=1.5)
	legend("bottomright",inset=c(0.01,0.01),legend=c("1 = dark selective","-1 = light selective"),xpd=TRUE,cex=1.5)

	selec_list[[1]]<-sapply(param_list, function(x) sum(x$counts[x$mid>=0.5]))
	selec_list[[2]]<-sapply(param_list, function(x) sum(x$counts[x$mid<=-0.5]))
	return(selec_list)


}

