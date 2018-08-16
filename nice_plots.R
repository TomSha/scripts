plot_clusters <- function(widths,file=NA,ext="ctrl_CD31"){
	mctrl=apply(data_all[,grep("ctrl",names(data_all))],1,mean)
	sdctrl=apply(data_all[,grep("ctrl",names(data_all))],1,sd)
			   
	if(!is.na(file)) pdf(file,10,5)
		par(mar=c(12,2,2,0))
		layout(matrix(1:2,ncol=2),widths=widths);
		clsub=cl[threshold]
		datasub=data[threshold,2:ncol(data)];
		#sets a side bar with each colour segment the length of each cluster 
		image(y=1:sum(threshold),z=matrix(clsub[order(clsub)],nrow=1),col=rainbow(max(cl)),xaxt='n',yaxt='n',ylab='n')
		#gets the centre of each cluster segment so that you can add labels (e.g. cluster 1)
		pos=floor(by(1:sum(threshold),clsub[order(clsub)],mean))
		axis(2,at=pos,labels=clsub[order(clsub)][pos],ylab=F)
		par(mar=c(12,0,2,2))
		rangecol=range(datasub);
		image(x=1:nEP,y=1:sum(threshold), z=t(as.matrix(datasub)[order(clsub),]),col=colorRampPalette(colors=c("blue","white","red"))(100),xaxt='n',yaxt='n',xlab="",ylab="",frame=F,breaks=seq(rangecol[1],rangecol[2],length.out=101))
		axis(1,at=1:nEP,labels=labels,las=2)
															      
		axis(1,at=1:length(external),labels=names(external),las=2)

		if(!is.na(file)) dev.off()
}
