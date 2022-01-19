
#################################################################
# Function:  Ramanome baseline function
# Author: Rongze Chen, Yuehui He
# Last update: 2021-11-01, Rongze Chen, Yuehui He
#################################################################
# install necessary libraries

## Clean R environment
rm(list=ls())
setwd("./")
sourcedir <- Sys.getenv("RamEx")
#sourcedir <- paste("/home/gene/jinggc/RamEX/databases/QC")
source(paste(sourcedir, "/databases/QC/Function/RamanEx_plot.R", sep = "/"))
source(paste(sourcedir, "/databases/QC/Function/folder_reader.R", sep = "/"))
source(paste(sourcedir, "/databases/QC/Function/spc_melt.R", sep = "/"))
source(paste(sourcedir, "/databases/QC/Function/util_ramex.R", sep = "/"))
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--raman_data_folder"),type="character", help="Input_the_raman_data_floder"),
	make_option(c("-o", "--out_dir"), type="character", default='Baseline_result', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$raman_data_folder)) stop('Please input a test raman file')

ramanpath <- opts$raman_data_folder 
outpath <- opts$out_dir


#outputpath creation
dir.create(outpath)
options(warn=-1)

Name_group <- c("group_A","group_B","group_C") #"C","D","E"   
all_spc <- folder_reader(ramanpath)
colnames(all_spc) <- gsub("X","",colnames(all_spc))
filename <- all_spc$filename
meta_data <-  as.data.frame(tstrsplit(gsub(".txt","",filename),"_"),col.names = Name_group)

all_data <- cbind(filename,meta_data,all_spc[,-1])
spectrum <- as.numeric(colnames(all_spc[,-1]))
data_hyperSpec<-new ("hyperSpec", data = data.frame(all_data[,1:7]),spc = all_data[,11: ncol(all_data)],wavelength=as.numeric(as.character(colnames(all_data[,11: ncol(all_data)]))))
data_hyperSpec
data_baseline_plot <- Spec_Nor_all(data_hyperSpec,c("group_A","group_B","group_C"))
fwrite(meta_data,paste(outpath,"meta_data.txt",sep = "/"),row.names = F,col.names = T,quote = F,sep = "\t")
ggsave(filename=paste(outpath,"/","baseline_nor-max-all-2.png", sep=""),plot=data_baseline_plot, limitsize=T,width = 16,height = 4)


