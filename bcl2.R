bcl <- function(y,wdt,lam,varB,varC,varX){
	N=length(y);
	B=rep(NA,N);
	c=rep(0,N);
	sks=rep(0,N);
	B[1]=mean(y);
	loglik=0;
	for(t in 2:N){
		cnew=c[t-1]*exp(-lam);
		Bnew=(varX*B[t-1]+varB*(y[t]-cnew))/(varX+varB);
		logp0=log(1-wdt)-1/2*log((2*pi)^2)-1/2*log(varX)-1/2*log(varB)-(Bnew-B[t-1])^2/(2*varB)*(1+varX/varB);
		cspike=(varC*(y[t]-B[t-1])+cnew*(varX+varB))/(varX+varC+varB); 
		Bspike=B[t-1]+varB/varC*(cspike-cnew);
		logp1=log(wdt)-1/2*log((2*pi)^3)-1/2*log(varX)-1/2*log(varB)-1/2*log(varC)-(Bspike-B[t-1])^2/(2*varB)*(1+varX/varB+varC/varB);
		
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
	}
	return(list(c,B,sks,loglik));
}

plot.tbcl <- function(y,nfdecay=10,cex=.1,pch=21,varB,varC,varX,pp){
	t<-y;
	dt<-bcl(t,pp,log(2)/nfdecay,varB,varC,varX);

	layout(matrix(1:2));
	par(mar=c(2,2,1,1));
	plot(t,type='l',xaxt='n',xlab="",ylab="", main=paste("cell",i));
	lines(dt[[2]],col="blue",lwd=3);
    
	plot(dt[[1]]/dt[[2]],type='l',xaxt='n',xlab="",ylab="");
	points(x=(1:length(t))[dt[[3]]==1],y=rep(0,length(t))[dt[[3]]==1],col="red",cex=cex,pch=pch);
}

bcl.EM<-function(y,varB0,varC0, varX0,pp0=0.4,lam0=log(2)/10){
	layout(matrix(1:2,ncol=1))
	pp=runif(n=1,min=0,max=pp0);
	varC=abs(rnorm(n=1,mean=varC0,sd=varC0));
	varB=abs(rnorm(n=1,mean=varB0,sd=varB0));
	varX=abs(rnorm(n=1,mean=varX0,sd=varX0));
	cat(varB,varC,varX,pp,'\n')
	lam=lam0
	par(ask=T)
	cat("Probe_decay"," ","noise"," ","Baseline_variance"," ","Calcium_variance"," ","spike_frequency"," ","log-likelihood","\n")
	oldvarX=varX;oldvarB=varB;oldvarC=varC;
	for(i in 1:100){
		trace=bcl(y,pp,lam,varB,varC,varX);
		cat(lam,' ',varX,' ',varB,' ',varC,' ',pp,' ',trace[[4]],"\n")
		cat(lam,' ',mean((y-trace[[1]]-trace[[2]])^2),' ',
		            mean(diff(trace[[2]])^2),' ',
					mean(((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam))[trace[[3]][2:length(y)]==1])^2),' ',pp,' ',trace[[4]],"\n");


		plot(y,type='l')
		lines(trace[[2]]+trace[[1]],col="red")
		lines(trace[[2]],col="blue")
		plot(trace[[1]],type='l')
	    points(x=(1:length(y))[trace[[3]]==1],y=rep(0,length(y))[trace[[3]]==1],col="red",cex=.2,pch=19);
		oldvarX=varX; oldvarB=varB; oldvarC=varC;
		varX=mean((y-trace[[1]]-trace[[2]])^2)
		varB=mean(diff(trace[[2]])^2)
		if(sum(trace[[3]])>0) varC=mean(((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam))[trace[[3]][2:length(y)]==1])^2)
		pp=sum(trace[[3]])/length(y)
	}
}

bcl.ml<-function(y,varB0,varC0, varX0,pp0=0.4,lam0=.001){
	#q90=quantile(y,0.5)
	#varX0=var(y[y<q90])
	layout(matrix(1:2,ncol=1))
	pp=pp0
	varC=varC0
	varB=varB0
	varX=varX0
	lam=lam0
	par(ask=T)
	cat("Probe_decay"," ","noise"," ","Baseline_variance"," ","Calcium_variance"," ","spike_frequency"," ","log-likelihood","\n")
	for(i in 1:30){
		trace=bcl(y,pp,lam,varB,varC,varX);
		cat(lam,' ',varX,' ',varB,' ',varC,' ',pp,' ',trace[[4]],"\n")
		plot(y,type='l')
		lines(trace[[2]]+trace[[1]],col="red")
		lines(trace[[2]],col="blue")
		plot(trace[[1]],type='l')
	    points(x=(1:length(y))[trace[[3]]==1],y=rep(0,length(y))[trace[[3]]==1],col="red",cex=.2,pch=19);
		
		#pp=sum(trace[[3]])/(length(y))
		varB_s1=mean(((trace[[2]][2:length(y)]-trace[[2]][1:(length(y)-1)])[trace[[3]][2:length(y)]==1])^2)
		varB_s0=mean(((trace[[2]][2:length(y)]-trace[[2]][1:(length(y)-1)])[trace[[3]][2:length(y)]==0])^2)
		varB_tot=mean(((trace[[2]][2:length(y)]-trace[[2]][1:(length(y)-1)]))^2)
		
		#varX=varX*exp(0.01*i)
		#varC=varC*exp(0.01*i)
		#varB=0.5*(varB_tot+sqrt(varB_tot^2+4*(2*(1-pp)*varX*varB_s0+2*pp*(varX+varC)*varB_s1)))
		if(sum(trace[[3]])>0) varC=mean(((trace[[1]][2:length(y)]-trace[[1]][1:(length(y)-1)]*exp(-lam))[trace[[3]][2:length(y)]==1])^2)
	}
}
