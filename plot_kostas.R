prefix0="2018-08-23-DotSpeed3Deg-BG-black_F1T3"
#prefix=paste(prefix0,"_rslice3",sep="")
folder=paste("~/Desktop/2018-08-23-DotSpeed3Deg-BG-black_F1T3/IMG/IMG_slice_1/results/2018-08-23-DotSpeed3Deg-BG-black_F1T3")
epfile=paste("~/Desktop/2018-08-23-Imaging_SpeedDot_F1_TimeRecords/2018-08-23-DotSpeed3Deg-BG-black_F1T3_time.log")
ep=read.table(epfile);
nEP<-nrow(ep)/2
ep$V1=floor((ep$V1+90)*20)
#source("bcl.R")

#ncells=system(paste("/home/meyer-lab/ToolBox_v3.6/bin/getcells ",folder,prefix,sep=""),intern=T);
read.cell <- function(id,xmax=700){
#	prefix<-paste(prefix0,"_rslice",slice,sep="")
	print(prefix)
#	layout(matrix(1:2))
    par(mar=c(10,3,2,1))
	cat(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder,sep=""))
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder,sep=""));
	t<-read.table("out_resp.dat")$V1;
	system(paste("/home/meyer-lab/ToolBox_v3.6/bin/detrend out_resp.dat out_dt.dat",length(t),"20 5"));
	t_dt<-read.table("out_dt.dat")$V1;
	#t_bcl<<-bcl(t,0.1,log(2)/30,1e-3,1);
	plot((t),type='l',ylim=(c(1,xmax)),main=paste("cell id",id),xaxt='n',ann=F);
	axis(1,at=ep$V1[1:nEP*2-1],label=ep$V2[1:nEP*2-1],las=2)
	rect(xleft=ep$V1[1:nEP *2 -1],xright=ep$V1[1:nEP *2 ],ybottom=rep(-1000,nEP*2),ytop=rep(10000,nEP*2),col=rgb(0,1,0,.1))
    par(mar=c(1,3,2,1))
#	plot(t_dt,type='l')
	#plot(t_bcl[[1]]/t_bcl[[2]],type='l');
#	points(x=(1:length(t))[t_bcl[[3]]==1],y=rep(0,length(t))[t_bcl[[3]]==1],col="red",cex=.2,pch=21)
}

