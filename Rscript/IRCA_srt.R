
#################################################################
# Function:  Ramanome IRCA spec range trimming function  
# Author: Yuehui He, Shi Huang
# Last update: 2021-09-22, Yuehui He,Shi Huang
#################################################################
# install necessary libraries

## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util_clean.R',sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--wnt_data"),type="character", help="Input_the_wave_number_trimming_data"),
	make_option(c("-o", "--out_dir"), type="character", default='spec_range_trimming', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$wnt_data)) stop('Please input the wave number trimming data')

matrixfile <- opts$wnt_data #"approxfun_Alldata.txt"
outpath <- opts$out_dir #"outputpath" 


#outputpath creation
dir.create(outpath)


options(warn=-1)
#-------------------------------
# Spec range trimming
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t");
mat <- data.frame(mat)
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
#-------------------------------
# scale by dividing sum area
#-------------------------------
rraw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-round(raw_wn, 0)
write.table(mat,paste(outpath,"/srt.table",sep=""),row.names=T,quote=F,sep="\t")
cat("The number of Raman shifts: ", ncol(mat) , "\n")

