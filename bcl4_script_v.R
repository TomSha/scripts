bcl<- function(y,wdt,lam,varB,varC,B0,Cmean,freq=5,mode="GAUS"){
	varX=var(diff(y)[diff(y)<0]);
	N=length(y);
	B=rep(NA,N);
	c=rep(0,N);
	sks=rep(0,N);
    B[1]=B0;
	loglik=0;
	dt=1./freq;
	for(t in 2:N){
		cnew=c[t-1]*exp(-lam*dt);
		Bnew=(varX*B[t-1]+varB*dt*(y[t]-cnew))/(varX+varB*dt);
		if(mode=="GAUS"){
			logp0=log(1-wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=Cmean+cnew+(y[t]-cnew-B[t-1])/(1+varB*dt/varC+varX/varC)
			Bspike=B[t-1]+varB*dt/varC*(cspike-cnew-Cmean);
			logp1=log(wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt+varC)-(y[t]-cnew-B[t-1]-Cmean)^2/(2*varX+2*varB*dt+2*varC);
		} else if (mode=="EXP"){
			logp0=log(1-wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=y[t]-B[t-1]-varC*(varB*dt+varX);
			Bspike=B[t-1]+varC*varB*dt;
			logp1=log(wdt)-0.5*log((2*pi)^2)-1/2*log(varX)-1/2*log(varB*dt)+log(varC)-(Bspike-B[t-1])^2/(2*varB*dt)*(1+varX/varB/dt) - varC*(cspike-cnew) ;

		} else if (mode=="EXP2"){
			prevs=sks[t-1]+1;
			logp0=log(1-wdt[prevs])-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=y[t]-B[t-1]-varC*(varB*dt+varX);
			Bspike=B[t-1]+varC*varB*dt;
			logp1=log(wdt[prevs])-0.5*log((2*pi)^2)-1/2*log(varX)-1/2*log(varB*dt)+log(varC)-(Bspike-B[t-1])^2/(2*varB*dt)*(1+varX/varB/dt) - varC*(cspike-cnew) ;
		}
		
		if(mode=="GAUS"){
			if(logp1<logp0){
				c[t]=cnew;
				B[t]=Bnew;
				loglik=loglik+logp0
			} else {
				c[t]=cspike;
				B[t]=Bspike;
				sks[t]=1;
				loglik=loglik+logp1
			}
		} else if(mode=="EXP" || mode=="EXP2"){
			if(logp1<logp0 || cspike-cnew<=0){
				c[t]=cnew;
				B[t]=Bnew;
				loglik=loglik+logp0;
			} else {
				c[t]=cspike;
				B[t]=Bspike;
				sks[t]=1;
				loglik=loglik+logp1
			}
		}
	}
	
	return(list(c,B,sks,loglik));
}

plot.tbcl <- function(y,pp,lambda=0.6,varB,varC,B0,Cmean,freq,mode="GAUS",lay=TRUE,cex=.1,pch=21){
	t<-y;
	dt<-bcl(y=t,
	        wdt=pp,
			lam=lambda,
			varB=varB,
			varC=varC,
			B0=B0,
			Cmean=Cmean,
			freq=freq,
			mode=mode);

    if(lay)	layout(matrix(1:2));
	par(mar=c(2,2,1,1));
	plot(t,type='l',xaxt='n',xlab="",ylab="");
	lines(dt[[2]],col="blue",lwd=3);
    
	plot(dt[[1]],type='l',xaxt='n',xlab="",ylab="");
	points(x=(1:length(t))[dt[[3]]==1],y=rep(0,length(t))[dt[[3]]==1],col="red",cex=cex,pch=pch);
}

bcl.EM<-function(y,pp0,lam0=0.6,varB0,varC0, B0=1 ,Cmean=0, freq=5,mode="GAUS"){

	zC=length(y)*0.2; zB=0.5; zX=1;
	layout(matrix(1:2,ncol=1))
	pp=pp0; #runif(n=1,min=0,max=pp0);
	
	varC=varC0;#abs(rnorm(n=1,mean=varC0,sd=varC0));
	varB=varB0;#abs(rnorm(n=1,mean=varB0,sd=varB0));

	cat(varB,varC,pp,'\n')
	lam=lam0;
	dt=1/freq;
	N=length(y)
	par(ask=T)
	cat("Probe_decay"," ","noise"," ","Baseline_variance"," ","Calcium_variance"," ","spike_frequency"," ","log-likelihood","\n")
	for(i in 1:100){
		trace=bcl(y,pp,lam,varB,varC,B0,Cmean,freq,mode);
		cat(lam,' ',varB,' ',varC,' ',pp,' ',Cmean,trace[[4]],"\n")
		cat(lam,' ',
		            mean(diff(trace[[2]])^2),' ',
					mean(((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam*dt)-Cmean)[trace[[3]][2:length(y)]==1])^2),' ',
					pp,' ',
					mean((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam*dt))[trace[[3]][2:length(y)]==1]),' ',
					trace[[4]],"\n");


		plot(y,type='l')
		lines(trace[[2]]+trace[[1]],col="red")
		lines(trace[[2]],col="blue")
		plot(trace[[1]],type='l')
	    points(x=(1:length(y))[trace[[3]]==1],y=rep(0,length(y))[trace[[3]]==1],col="red",cex=.2,pch=19);
		cat("Calculations\n")

		cat("varX=",varX,"\n");
		varB=(sqrt(N^2+8*sum(diff(trace[[2]])^2/dt))-N)/(4 *zB)
		cat("varB=",varB,"\n");
		B0=mean(trace[[2]])
		if(sum(trace[[3]])>0){
			if(mode=="GAUS"){
				Cmean=mean((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam*dt))[trace[[3]][2:length(y)]==1]);
				varC=mean(((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam*dt)-Cmean)[trace[[3]][2:length(y)]==1])^2)
			} else if(mode=="EXP") {
				varC=sum(trace[[3]])/(zC+sum(trace[[3]][2:N]
											 *(trace[[1]][2:N]-trace[[1]][1:(N-1)]*exp(-lam*dt))))
				cat("gammaC=",varC,"\n");
			}
		pp=sum(trace[[3]])/N
		} else {
			pp=1./N;
		}
	}
}


#generates a correlation value for each cell, where 1= always firing during epoch and silent outside of epochs, -1=always firing outside of epochs and never firing within epochs
#args <- commandArgs(trailingOnly= TRUE)


EpochCorr<-function(prefix,resultsfolder){ 
#added in a random noise model to calculate threshold for each cell
#thresh = 0.625 for dots
#thresh = 0.44 for simple
	
	

	bin="/home/meyer-lab/ToolBox_v3.6/bin/"
	binary<-paste(resultsfolder,"/",prefix,".bin",sep="");
	
	#if analysing volumes removes the slice # from prefix
	if (grepl("_rslice",prefix)){
		gprefix<-gsub("_rslice.*",x=prefix,replacement="")
		epochfile<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",gprefix,"_time.log",sep="");
		print(epochfile)
	}else{
	epochfile<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="");
	}
	outputC<-read.table(paste(resultsfolder,"/",prefix,"_outputC.dat",sep=""));
	outputRC<-read.table(paste(resultsfolder,"/",prefix,"_output_rawC.dat",sep=""));
	
#	thresh<-as.numeric(thresh)
#	print(thresh)

	system(paste(bin,"/readall ",binary,sep=""));	
	allresp<-read.table("out_resp_all.dat");
	nz=ncol(allresp);
	ncells=nrow(allresp);
	ep<-read.table(epochfile);
	matep<-matrix(floor(4*(ep$V1+90)),ncol=2,byrow=T);
	nEP=nrow(ep)/2;
	corvec<-rep(0,ncells);
	thresh<-rep(0,ncells);
	stim<-rep(-1,nz);
	for(i in 1:nEP) stim[matep[i,1]:matep[i,2]]=1;

	for(i in 1:ncells){
		trace<-as.numeric(allresp[i,])
		bb<-bcl(trace/mean(trace),0.01,0.8,9e-3,1e-3,.4,0,20,"EXP");
		ss<-bb[[3]]; ss[ss==0]=-1;
		corvec[i]=sum(ss*stim)/nz

		#calculate threshold with random noise model
#		p<-table(ss)[2]/table(ss)[1]
		p<-table(ss)[2]/length(ss)
		#mean = frames(epochs)-frames(outepochs)/frames * (2p-1)
		threshMean<-(sum(stim)/length(stim))*(2*p-1)

		#variance = 1/N 4p(1-p)
		threshSD<-sqrt(1/length(stim)*4*p*(1-p))

		thresh[i]<-threshMean+(4*threshSD)

	}


outputCT<-outputC[corvec>thresh,];
outputRCT<-outputRC[corvec>thresh,];


	print(paste("outputC" , nrow(outputC), sep=""))
	print(paste("outputCT ", nrow(outputCT), sep=""))

write.table(corvec,paste(resultsfolder,"/",prefix,"_corvec.dat",sep=""),row.names=F,col.names=F);
write.table(outputCT,paste(resultsfolder,"/",prefix,"_outputCT.dat",sep=""),row.names=F,col.names=F);
write.table(outputRCT,paste(resultsfolder,"/",prefix,"_output_rawCT.dat",sep=""),row.names=F,col.names=F);
write.table(corvec>thresh,paste(resultsfolder,"/",prefix,"_noise_thresh.dat",sep=""),row.names=F,col.names=F);

}


EpochCorr2 <- function(prefix){
		
		resultsfolder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/im/im_slice_1/results_toolbox/",sep="")

		bin="/home/meyer-lab/ToolBox_v3.6/bin/"
		binary<-paste(resultsfolder,prefix,".bin",sep="")
		epochfile<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="")

	    system(paste(bin,"/readall ",binary,sep=""));
	    allresp<-read.table("out_resp_all.dat");
	    nz=ncol(allresp);
	    ncells=nrow(allresp);
	    ep<-read.table(epochfile);

		#take from 1 to 1 before the start frame
		start<-floor(20*(ep[1,1]+90))-1
	    nEP=nrow(ep)/2;
	    corvec<-rep(0,ncells);

		for(i in 1:ncells){
			trace<-as.numeric(allresp[i,])
			bb<-bcl(trace/mean(trace),0.01,0.8,9e-3,1e-3,.4,0,20,"EXP");
			ss<-bb[[3]][1:start]
			ss[ss==0]=-1;
			corvec[i]=sum(ss*-1)/start
			
			}
			
		write.table(corvec,paste(resultsfolder,"/",prefix,"_corvecPreB.dat",sep=""),row.names=F,col.names=F);

}


plot_cell<-function(id,prefix0,slice){
	
	prefix<-paste(prefix0,"_rslice",slice,sep="")
	folder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix0,"/results_toolbox/",sep="");
	epfile<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix0,"_time.log",sep="")
	ep<-read.table(epfile)$V1
	ep<-floor((ep+97)*4)
	nEP<-length(ep)/2
	bin="/home/meyer-lab/ToolBox_v3.6/bin/"

	print(paste(bin,"/read"," ",id," ",folder,prefix,"\n",sep=""));
	system(paste(bin,"/read"," ",id," ",folder,prefix,"\n",sep=""));
	t<-read.table("out_resp.dat")$V1;
#	t<-bcl(t/mean(t),0.01,0.2,1e-3,1e-3,.4,0,20,"EXP");
	t<-bcl(t/mean(t),0.01,0.8,9e-3,1e-3,.4,0,20,"EXP");
	spikes<-t[[3]];
	nz<-length(t[[3]]);
	spikes<-which(spikes==1);
	yax<-rep(-3,length(nz));
	spikes<-cbind(spikes,yax);
#	par(mfrow=c(2,1));
#	par(mar=c(0,3,2,2))
#	plot(t[[2]],type="l",ann=F,yaxt="n",xaxt="n",bty="n",ylim=c(-4,max(t[[2]])));
	par(mar=c(0,3,0,2))
	plot(t[[1]],type="l",ann=F,yaxt="n",xaxt="n",bty="n",ylim=c(-4,max(t[[1]])));
	points(spikes,pch=20);
	rect(xleft=ep[1:nEP*2-1],xright=ep[1:nEP*2],ybottom=rep(0,nEP*2),ytop=rep(max(t[[1]]),nEP*2),col=rgb(0,1,0,.1))


}


plot_hist<-function(vec){

	if (min(vec)<0){
		xlims<-c(-0.1,0.4)
		} else {
		xlims<-c(0.45,0.75)
		}
	
	minbreak<-min(vec)
	maxbreak<-max(vec)

	hist(vec,freq=F,xlim=xlims,col=rgb(0,0,1,0.3),breaks=seq(minbreak,maxbreak,l=20))

	}



args<-commandArgs(trailingOnly = T);
EpochCorr(args[1],args[2]);

