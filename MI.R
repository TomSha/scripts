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
infotab=rep(0,ncells);


for (cell in 1:ncells){

	system(paste(bin,"/read ",cell," ",folder,"/",prefix,sep=""));
	trace<-read.table("out_resp.dat")$V1;

	f=min(trace);
	t=max(trace);

    Dist<-density(x=trace[MatEP[1,1]:MatEP[nEP,2]], from=f,to=t,n=ngrid);
	Dist$y<-Dist$y/sum(Dist$y);

	NullDist<-density(x=trace[1:MatEP[1,1]], from=f,to=t,n=ngrid)
	NullDist$y<-NullDist$y/sum(NullDist$y);
    
	Info=0;
	for(j in 1:ngrid){
		if(NullDist$y[j]>1e-10){
			Info<-Info+NullDist$y[j]*log(NullDist$y[j]/Dist$y[j])/log(2.);
		}
	}

	infotab[cell]=Info;

}

write.table(infotab,file="KL.dat",row.names=F,col.names=F)

