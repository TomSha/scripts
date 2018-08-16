require(plotrix)
require(mmand)

pulledresultsfolder=(paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/pulled_results_volumes/",sep=""));
pulledoutput=read.table(paste(pulledresultsfolder,"/","all_centers.dat",sep=""));
pulledinc=read.table(paste(pulledresultsfolder,"/","pulled_inc.dat",sep=""))
pulleddencl=read.table(paste(pulledresultsfolder,"/","pulled_summary.dat",sep=""))
pulleddencl_sort=pulleddencl[order(pulleddencl$V1),];
#pulleddencl_inc=pulleddencl[pulledinc[,1],];
pulledoutput_inc=pulledoutput[pulledinc[,1],];
toolbox=paste("/home/meyer-lab/ToolBox_v3.6/")
prefix_list=read.table(paste(toolbox,"/scripts/prefix_list_vol",sep=""))$V1

epochsC_list<-vector("list",length(prefix_list));
output0_list<-vector("list",length(prefix_list));
outputR0_list<-vector("list",length(prefix_list));
dencl_list<-vector("list",length(prefix_list));
inc_list<-vector("list",length(prefix_list));
output_list<-vector("list",length(prefix_list));
outputR_list<-vector("list",length(prefix_list));
epochsC_list<-vector("list",length(prefix_list));
cellcen_list<-vector("list",length(prefix_list));

for (i in 1:length(prefix_list)){
	prefix<-prefix_list[i];
	resultsfolder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/H2Bslow/",prefix,"/results_toolbox",sep="");
    timelogfolder<-paste("/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/",prefix,"_time.log",sep="");

	epochsC_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_rslice10_epochsC.dat",sep=""));
	ord<-order(epochsC_list[[i]]);
	epochsC_list[[i]]<-epochsC_list[[i]][,1][ord]

	output0_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_combined_outputC.dat",sep=""));
	outputR0_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_combined_output_rawC.dat",sep=""));
	output0_list[[i]][,3:ncol(output0_list[[i]])]<-output0_list[[i]][,3:ncol(output0_list[[i]])][ord];
	outputR0_list[[i]][,3:ncol(outputR0_list[[i]])]<-outputR0_list[[i]][,3:ncol(outputR0_list[[i]])][ord];


	dencl_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_summary.dat",sep=""));
	dencl_list[[i]]<-dencl_list[[i]][order(dencl_list[[i]]$V1),];
	inc_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_inc.dat",sep=""));
	dencl_list[[i]]<-dencl_list[[i]][inc_list[[i]][,1],];
	
	output_list[[i]]<-output0_list[[i]][output0_list[[i]]$V2 %in% dencl_list[[i]]$V1,];
	outputR_list[[i]]<-outputR0_list[[i]][outputR0_list[[i]]$V2 %in% dencl_list[[i]]$V1,];
	output_list[[i]]<-setNames(output_list[[i]],c("V1","V2",as.character(epochsC_list[[i]])));

	cellcen_list[[i]]<-read.table(paste(resultsfolder,"/",prefix,"_cell_centers.dat",sep=""));
}


normoutput<-vector("list",length(prefix_list));
selec<-vector("list",length(prefix_list));
ent<-vector("list",length(prefix_list));

maxn<-function(x) sort(x,decreasing = T)[1];

for (k in 1:length(prefix_list)){

	#calucluate max epoch value
	maxepochvalue<-apply(output_list[[k]][,3:ncol(output_list[[k]])],1,maxn);
	maxepochname<-colnames(output_list[[k]])[max.col(output_list[[k]][,3:ncol(output_list[[k]])])+2]
	#calculate entropy
	normoutput[[k]]<-outputR_list[[k]][,3:ncol(outputR_list[[k]])]/rowSums(outputR_list[[k]][,3:ncol(outputR_list[[k]])])
	ent[[k]]<-apply(normoutput[[k]],1,function(x) sum(x*log2(1/x)))
	
    no<-nrow(output_list[[k]])
    selec[[k]]<-data.frame(maxepochname=rep(NA,no),entropy=rep(NA,no),slice=rep(NA,no),cellID=rep(NA,no),count=rep(NA,no),C1=rep(NA,no),C2=rep(NA,no))
    for (m in 1:nrow(output_list[[k]]))
    {
		prefix=prefix_list[k]
		selec[[k]]$maxepochname[m]  <- maxepochname[m]
		selec[[k]]$entropy[m]       <- ent[[k]][m];
		selec[[k]]$slice[m]			<- output_list[[k]][m,1];
		selec[[k]]$cellID[m]        <- output_list[[k]][m,2]-sum(cellcen_list[[k]]$V3<output_list[[k]][m,1]);
		selec[[k]]$count[m]			<- output_list[[k]][m,2]
		selec[[k]]$C1[m]            <- dencl_list[[k]]$V6[m]+1
		if(any(pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m])){
			selec[[k]]$C2[m] <- pulleddencl_sort$V6[pulledinc[,1] & pulledoutput$V2==prefix & pulledoutput$V4==selec[[k]]$C1[m]]+1
		} else {
	        selec[[k]]$C2[m] <-"NULL"
																																	  			}
		}			
	}

##FUNCTIONS##


#####print cells of a specific cluster

cells_in_cluster <-function(cluster){
	a <- vector("list",length(prefix_list))
	percent<-   "NULL"
	cell_num <- "NULL"
	for (j in 1:length(prefix_list)){
		a[[j]]<-selec[[j]][which(selec[[j]][,7]==cluster),]
		cell_num[j]<-nrow(a[[j]])                       }
  		print(a)    
		print(cell_num)
		sum(as.integer(cell_num))
		}
#####
