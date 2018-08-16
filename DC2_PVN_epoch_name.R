args <- commandArgs(trailingOnly = TRUE);

epochs_list<-args[1]
epochs_list <- read.table(epochs_list)
epochs_list<-apply(epochs_list, 1, function(x) paste(x,"_time.log",sep=""))

#epochs_list <- read.table("/home/meyer-lab/SeG/R/epoch_list.txt")
epochs_folder <- "/media/meyer-lab/Elements/Work/Stimulus_Barrage/barrage_v2/ExpLogs_TomSh/"

for (i in 1:length(epochs_list))
{
	epochs <- read.table(paste(epochs_folder,epochs_list[i],sep=""))

	epochs_new <- gsub("LEFT","90",epochs$V2)
	epochs_new <- gsub("RIGHT","270",epochs_new)
	epochs_new <- gsub("UP","0",epochs_new)
	epochs_new <- gsub("DOWN","180",epochs_new)
	epochs_new <- gsub("R_L","R_90",epochs_new)
	epochs_new <- gsub("R_R","R_270",epochs_new)
	epochs_new <- gsub("R_U","R_0",epochs_new)
	epochs_new <- gsub("R_D","R_180",epochs_new)
	epochs_update <-cbind(epochs$V1,epochs_new)

	#write(epochs,(paste(epochs_folder,epochs_list[i,],sep="")))
	file <- paste(epochs_folder, epochs_list[i], sep="")
	write.table(epochs_update,file, row.names = F, col.names=F, quote=F)

	}



