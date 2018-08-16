args <- commandArgs(trailingOnly=T);

bin=args[1]
folder=args[2]
filename=args[3]
ep=read.table(args[4])
MatEP<-matrix(floor((ep$V1+90)*20),ncol=2,byrow=T);
nEP=nrow(ep)/2;
ngrid=100;
prefix=strsplit(filename,split="[.]")[[1]][1]

print(paste(bin,"/getcells ",folder,"/",filename,sep=""));
ncells=as.numeric(system(paste(bin,"/getcells ",folder,"/",filename,sep=""),intern=T));
cat(nEP,ncells,"\n");
infotab=matrix(0,ncol=nEP,nrow=ncells);


for (cell in 1:ncells){

	system(paste(bin,"/read ",cell," ",folder,"/",prefix,sep=""));
	trace<-read.table("out_resp.dat")$V1;

	f=min(trace);
	t=max(trace);

	dlist<-list();

	for(i in 1:nEP){
		dlist[[i]]<-density(x=trace[MatEP[i,1]:MatEP[i,2]],
		                    from=f,to=t,n=ngrid);
		dlist[[i]]$y<-dlist[[i]]$y/sum(dlist[[i]]$y);
	}

	NullDist<-density(x=trace[1:MatEP[1,1]], from=f,to=t,n=ngrid)
	NullDist$y<-NullDist$y/sum(NullDist$y);

	TotDist<-rep(0,ngrid);
	for(i in 1:ngrid){
		for(k in 1:nEP){
			TotDist[i]<-TotDist[i]+dlist[[k]]$y[i]/nEP;
		}
	}


	InfoVec<-rep(0,nEP);
	hist<-matrix(nrow=nEP,ncol=ngrid);

	for(i in 1:nEP){
		for(j in 1:ngrid){
			if(dlist[[i]]$y[j]>1e-10){
				InfoVec[i]<-InfoVec[i]+dlist[[i]]$y[j]*log(dlist[[i]]$y[j]/TotDist[j])/log(2.);
			}
		}
	}

	infotab[cell,]=InfoVec;

}

write.table(infotab,file="infotab.dat",row.names=F,col.names=F)

