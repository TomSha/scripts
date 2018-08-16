library(proxy)

prefix_list<-c("170203_PVN_barrage_7dpf_F1_1",
               "170203_PVN_barrage_7dpf_F1_2",
			   "170203_PVN_barrage_7dpf_F2_1",
			   "170203_PVN_barrage_7dpf_F2_2",
			   "170210_PVN_barrage_7dpf_F1_1",
			   "170210_PVN_barrage_7dpf_F1_2",
			   "170210_PVN_barrage_7dpf_F2_1",
			   "170210_PVN_barrage_7dpf_F2_2");

epochs_ord<-c("BAR_0","BAR_180","BAR_270","BAR_90","GRAT_0_SF15","GRAT_0_SF3","GRAT_0_SF30","GRAT_180_SF15","GRAT_180_SF3","GRAT_180_SF30","GRAT_270_SF15","GRAT_270_SF3","GRAT_270_SF30","GRAT_90_SF15","GRAT_90_SF3","GRAT_90_SF30","LIGHTBAR_0","LIGHTBAR_180","LIGHTBAR_270","LIGHTBAR_90","LOOM","WIGGLY_15","WIGGLY_3","WIGGLY_30");
data_list=vector("list",length(prefix_list));
exp_ID_list<-vector("list",length(prefix_list));
cen<-vector("list",length(prefix_list))
mip<-vector("list",length(prefix_list))
dims<-vector("list",length(prefix_list))
for(i in 1:length(prefix_list)) {
	data_list[[i]]<-as.matrix(read.table(paste("data/",prefix_list[i],"_output.dat",sep="")))[,2:25]
	elist<-read.table(paste("data/",prefix_list[i],"_time.log",sep=""))$V2
	elist<-elist[1:(length(elist)/2) *2 ]
	elist<-elist[grep("WIGGLY",elist,invert=T)];
	elist<-elist[grep("CONCENTRIC",elist,invert=T)];
	elist<-c(as.character(elist),"WIGGLY_3","WIGGLY_15","WIGGLY_30");
	data_list[[i]]<-data_list[[i]][,order(elist)];
	exp_ID_list[[i]]<-rep(prefix_list[i],nrow(data_list[[i]]))
	cen[[i]]<-read.table(paste("data/",prefix_list[i],"_cell_centers.dat",sep=""));
	mip[[i]]<-read.table(paste("data/",prefix_list[i],"_mip.dat",sep=""));
	dims[[i]]<-as.numeric(read.table(paste("data/",prefix_list[i],"_sizes.dat",sep="")));
}
data<-as.matrix(do.call(rbind,data_list))
colnorms<-apply(data,2,mean);
for(i in 1:ncol(data)) data[,i]=data[,i]/colnorms[i];
rownorms<-apply(data,1,mean);
for(i in 1:nrow(data)) data[i,]=data[i,]/max(data[i,]);
data<-t(apply(data,1,function(x) (x-min(x))/(max(x)-min(x))))

exp_ID<-unlist(exp_ID_list)
cen<-as.matrix(do.call(rbind,cen))

N=50

get_CH <- function(K,clcl){
	cutcl=cutree(clcl,k=K);
	global_center<-apply(data,2,mean);
	cl_centers<-matrix(NA,ncol=ncol(data),nrow=K);
	ns<-as.numeric(table(as.numeric(cutcl)));
	B=0; W=0;
	for(i in 1:K) {
		if(ns[i]==1) {
			cl_centers[i,]=data[which(as.numeric(cutcl)==i),];
		} else { 
			cl_centers[i,]=apply(data[which(as.numeric(cutcl)==i),],2,mean);
			W=W+sum((t(data[which(as.numeric(cutcl)==i),])-cl_centers[i,])^2);
		}
		B=B+ns[i]*sum((cl_centers[i,]-global_center)^2);
	}
	#return(c(B,W))
	return(W)
}

Get_sol<- function(shuf=F,DATA,N=50){
	if(shuf) {
		data=t(apply(DATA,1,sample));
	} else {
		data=DATA
	}
    distmat<-dist(data)
	clcl <- hclust(distmat)
	Wlist<- sapply(1:N,get_CH,clcl=clcl)
	mulist=seq(0,1000,0.1);
	tmp<-vector("numeric",length=N);
	for(i in 1:length(mulist)) tmp[i]=which.min(Wlist+mulist[i]*(1:N))
	return(tmp)
}

FreeEnergy <- function(mu){
	tmp<-vector("numeric",length=N);
	for(i in 1:N) tmp[i]=Wlist[i]+mu*i
    return(min(tmp));
}

FreeEnergy.K <- function(mu){
	tmp<-vector("numeric",length=N);
	for(i in 1:N) tmp[i]=Wlist[i]+mu*i
    return(which.min(tmp));
}

get_r <- function(n,N){
	r<-vector("list",n);
	r[[1]]=Get_sol(F,data,N);
	for(i in 2:(n+1)){
		r[[i]]=Get_sol(T,data,N)
		flush.console();
	    cat(paste("\r",i-1));
	}
	cat("\n");
	return(r)
}

rstat<-function(r,N){
	nsampling=length(r)-1
	lab<-matrix(0,ncol=N,nrow=nsampling);
	counter<-rep(0,N);
	av<-rep(NA,N);
	ss<-rep(NA,N);
	qq<-rep(NA,N);
	for(i in 1:nsampling){
		tab=table(r[[i+1]]);
		counter[as.numeric(names(tab))]=counter[as.numeric(names(tab))]+1;
		lab[i,as.numeric(names(tab))]=lab[i,as.numeric(names(tab))]+as.numeric(tab);
	}

	for(i in 1:N){
		tmp<-lab[,i]
		av[i]=mean(tmp[tmp>0]);
		ss[i]=sd(tmp[tmp>0]);
		qq[i]=as.numeric(quantile(tmp[tmp>0],prob=.9))
	}

	if(length(counter)==N) {
		return(list(av,ss,qq))
	} else { 
		return(NA);
	}
}

do_analysis1 <- function(n){
	r<<-get_r(n,50)
	rs<<-rstat(r,50)
	plot(table(r[[1]]))
	lines(rs[[1]]+rs[[2]],col="red")
}

do_analysis2 <- function(n){
	centers<<-matrix(NA,ncol=24,nrow=n);
	centers_sd<<-matrix(NA,ncol=24,nrow=n);
	distmat<-dist(data)
	hCL <<- hclust(distmat)
	ct<<-cutree(hCL,k=n)
	ordcl<<-rev(order(as.numeric(table(ct))));
	for(i in 1:n){
		if(is.null(nrow(data[ct==i,]))){
			centers[i,]<<-data[ct==i,];
			centers_sd[i,]=rep(0,24);
		} else {
			centers[i,]<<-apply(data[ct==i,],2,mean)
			centers_sd[i,]<<-apply(data[ct==i,],2,sd)
		}
	}	
}

show <- function(CLID,c="red",pch=19,cex=.6){
	layout(matrix(1:9,ncol=3));
	for(i in 1:8){
		image(1:dims[[i]][1],1:dims[[i]][2],matrix(mip[[i]]$V2,ncol=dims[[i]][2]),col=grey.colors(100))
		points(cen[ct==CLID & exp_ID==prefix_list[i],c(2,1)]+c(1.5,1.5),col=c,cex=cex,pch=pch)
	}
	plot_traj(CLID)
}

plot_traj<-function(X,min=0,max=1,NS=1,AS=1){
	nEP=24;
	par(las=2,mar=c(3,20,2,3))
	ord2<-order(centers[X,]);
	pos<-barplot(centers[X,ord2],names.arg=epochs_ord[ord2],horiz=T,xlim=c(0,centers[X,ord2[nEP]]+centers_sd[X,ord2[nEP]]),cex.axis=AS,cex.names=NS);
	arrows(y0=pos,y1=pos,x0=centers[X,ord2]-centers_sd[X,ord2],x1=centers[X,ord2]+centers_sd[X,ord2],code=3,length=0);
}

