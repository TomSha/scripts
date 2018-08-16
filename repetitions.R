reps1<-matrix(nrow=3,ncol=100)
reps2<-matrix(nrow=3,ncol=100)
reps3<-matrix(nrow=3,ncol=100)

reps1[1,]<-dnorm((1:100),mean=50,sd=4)
reps1[2,]<-dnorm((1:100),mean=25,sd=4)
reps1[3,]<-dnorm((1:100),mean=50,sd=8)

reps2[1,]<-dnorm((1:100),mean=50,sd=4)
reps2[2,]<-dnorm((1:100),mean=45,sd=4)
reps2[3,]<-dnorm((1:100),mean=50,sd=6)

reps3[1,]<-dnorm((1:100),mean=50,sd=4)
reps3[2,]<-dnorm((1:100),mean=50,sd=4)
reps3[3,]<-dnorm((1:100),mean=50,sd=4)

RQI<-function(X){
	x<-var(colMeans(X))/mean(apply(X,1,var))
	return(x)
	}

plot_reps<-function(X){
	par(mfrow=c(3,1));
	for (i in 1:nrow(X)) plot(X[i,],ylim=c(0,0.1),type="l",xlab="frame number",ylab="Fluorescence",main=paste("repetition ",i,sep=""));
}

plot_mean<-function(X){
	labs<-c("rep1", "rep2","rep3")
	BC<-barplot(rowMeans(X),xlab="mean response",ylab="Fluorescence")
	axis(1,labels=labs,at=BC)
}





#bin<-"/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/170203_PVN_barrage_7dpf_F1_1/im/im_slice_1/results_toolbox/170203_PVN_barrage_7dpf_F1_1.bin"
id<-1

#cat(paste("/home/meyer-lab/ToolBox_v3.6/bin/read"," ",id," ",folder,prefix,"\n",sep="")

