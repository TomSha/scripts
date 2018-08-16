bcl <- function(y,wdt,lam,rb,rs){
	N=length(y);
	B=rep(NA,N);
	c=rep(0,N);
	sks=rep(0,N);
	B[1]=mean(y);
	q90=quantile(y,0.9)
	varF=var(y[y<q90])
	for(t in 2:N){
		cnew=c[t-1]*exp(-lam);
		Bnew=(B[t-1]+rb*(y[t]-cnew))/(1+rb);
		p0=(1-wdt)*exp(-(y[t]-cnew-Bnew)^2/(2*varF)-(Bnew-B[t-1])^2/(2*rb*varF));
		#cspike=(rs*(y[t]-B[t-1])+cnew*(1-rs-rb))/(1+rs+rb); # This seems to be wrong!! 
		cspike=(rs*(y[t]-B[t-1])+cnew*(1+rb))/(1+rs+rb); 
		Bspike=B[t-1]+rb/rs*(cspike-cnew);
		p1=wdt*exp(-(y[t]-cspike-cnew-Bspike)^2/(2*varF)-(Bspike-B[t-1])^2/(2*rb*varF)-(cspike-cnew)^2/(2*rs*varF));
		
		if(p1<p0){
			c[t]=cnew;
			B[t]=Bnew;
		} else {
			c[t]=cspike+cnew;
			B[t]=Bspike;
			sks[t]=1;
		}
	}
	return(list(c,B,sks));
}

plot.bcl <- function(dataset,i,nfdecay=10,cex=.1,pch=21){
	t<-as.numeric(dataset[i,]);
	dt<-bcl(t,0.1,log(2)/nfdecay,1e-3,.1);
#dt<-bcl(t,0.001,log(2)/10,0.01,0.01,1)
	layout(matrix(1:2));
	par(mar=c(2,2,1,1));
	plot(t,type='l',xaxt='n',xlab="",ylab="", main=paste("cell",i));
	lines(dt[[2]],col="blue",lwd=3);
    
	plot(dt[[1]],type='l',xaxt='n',xlab="",ylab="");
	points(x=(1:length(t))[dt[[3]]==1],y=rep(0,length(t))[dt[[3]]==1],col="red",cex=cex,pch=pch);
}

plot.tbcl <- function(y,nfdecay=10,cex=.1,pch=21){
	t<-y;
	dt<-bcl(t,0.1,log(2)/nfdecay,1e-3,1);
#dt<-bcl(t,0.001,log(2)/10,0.01,0.01,1)
	layout(matrix(1:2));
	par(mar=c(2,2,1,1));
	plot(t,type='l',xaxt='n',xlab="",ylab="", main=paste("cell",i));
	lines(dt[[2]],col="blue",lwd=3);
    
	plot(dt[[1]]/dt[[2]],type='l',xaxt='n',xlab="",ylab="");
	points(x=(1:length(t))[dt[[3]]==1],y=rep(0,length(t))[dt[[3]]==1],col="red",cex=cex,pch=pch);
}
