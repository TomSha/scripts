#folder<-"/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/pulled_results_dots" #CHANGE HERE

#centers<-read.table(paste(folder,"/all_centers.dat",sep=""));
#centers_sd<-read.table(paste(folder,"/all_centers_sd.dat",sep=""));

prefix<-"180330_PVN_barrage_F1_1";
folder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"im/im_slice_1/results_toolbox",sep="");

centers<-read.table(paste(folder,prefix,"_outputCT.dat"))

cols<-1:(ncol(centers)-4);
numcol<-ncol(centers)-4;
numrow<-nrow(centers);

pcenters<-as.data.frame(matrix(0,ncol=ncol(centers),nrow=nrow(centers)));
#pcenters_sd<-as.data.frame(matrix(0,ncol=ncol(centers),nrow=nrow(centers)));


mat<-matrix(rep(cols,numrow),ncol=numcol,byrow=T);
pmat<-t(apply(mat,1,sample));

for (i in 1:nrow(pmat)){
	pcenters[i,5:ncol(pcenters)]<-centers[i,5:ncol(centers)][pmat[i,]]
#	pcenters_sd[i,5:ncol(pcenters_sd)]<-centers_sd[i,5:ncol(centers_sd)][pmat[i,]]
}
	
pcenters[,1:4]<-centers[,1:4];
#pcenters_sd[,1:4]<-centers_sd[,1:4];



#write.table(pcenters,file=paste(folder,"_permute/all_centers.dat",sep=""),row.names=F,col.names=F);
#write.table(pcenters_sd,file=paste(folder,"_permute/all_centers_sd.dat",sep=""),row.names=F,col.names=F);

write.table(pcenters,file=paste(folder,"_permute/outputCT.dat",sep=""),row.names=F,col.names=F);

