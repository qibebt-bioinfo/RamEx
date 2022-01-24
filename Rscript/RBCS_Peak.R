# Function:  Ramanome RBCS Peak function  
# Author: Shi Huang, Rongze Chen
# Last update: 2021-12-07, Rongze Chen, Gongchao Jing
#################################################################
# install necessary libraries
## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
sourcedir <- paste(sourcedir, "/databases/QC/Function",sep="")
source(paste(sourcedir, "util_ramex.R", sep = "/"))
source(paste(sourcedir, "peak_util.R", sep = "/"))
source(paste(sourcedir, "folder_reader.R", sep = "/"))
source(paste(sourcedir, "new_hyperSpec.R", sep = "/"))
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
        make_option(c("-i", "--raman_data_folder"),type="character", help="Input_the_raman_data_floder"),
	make_option(c("-m", "--meta_data"),type="character", help="Input_the_raman_data_RBCS_group_file"),
	make_option(c("-w", "--wave_num"),type="character", help="Input_the_raman_data_wave_number"),
        make_option(c("-f", "--data_format"), type="character", default="Horiba", help="Input_the_data_format[default %default]"),
        make_option(c("-p", "--perplexity"), type="character", default="5", help="Input_the_perplexity_of_the_tSNE[default %default]"),
        make_option(c("-o", "--out_dir"), type="character", default='Baseline_result', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$raman_data_folder)) stop('Please input a test raman file')
ramanpath <- opts$raman_data_folder
outpath <- opts$out_dir
data_format <- opts$data_format
meta_data <- opts$meta_data
perplexity <- as.numeric(opts$perplexity)
wave_num <- opts$wave_num
#outputpath creation
dir.create(outpath)


#读取数据并根据命名规则生成metadata文件 
Name_group <- c("group_A","group_B","group_C", "RBCS_group")

###数据标准化
if(data_format != "Horiba"){
        data_group <- renishaw(ramanpath,Name_group,save = T,outpath) # renishaw的数据
        fwrite(data_group,paste(outpath,"wavenumber_group.txt",sep = "/"),row.names = F,col.names = T,quote = F,sep = "\t")
        all_spc <- folder_reader(paste(outpath,"all_single_cell",sep = "/")) } else{
        all_spc <- folder_reader(ramanpath)
}
Meta_data <- read.table(meta_data, sep="\t",header=TRUE)
wave_num <- read.table(wave_num, sep="\t")


###平滑
#all_data <- S_G_smooth(all_spc)
data_hyperSpec <- new_hyperSpec(all_spc,Name_group)
data_hyperSpec$RBCS_group <- Meta_data$RBCS_group
#colnames(all_spc)
###质量控制
#data_hyperSpec_good <- spc_quanlity_control(data_hyperSpec,abs_noise = 0.1,sd_noise = 0.03)
#data_hyperSpec_good_plot <- Spec_all(data_hyperSpec_good)
#ggsave(filename=paste(outpath,"/","Original_all.png", sep=""),plot=data_hyperSpec_good_plot, limitsize=T,width = 16,height = 4)

###baseline

#data_baseline <- baseline(data_hyperSpec,poly.order = 7)
############
##baseline##
############
data_baseline <- data_hyperSpec-spc.fit.poly.below(data_hyperSpec, data_hyperSpec, poly.order = 7)
#plot (data_baseline, col = 1 : 2)
#write.csv(data_baseline,paste(output,"Cells_bg_baseline.csv",sep=""),quote = F,row.names = F)

########
##zero##
########
data_baseline_zero_hyperSpec <- data_baseline
#yminset <- apply (data_baseline$spc, 1, min)
#data_baseline_zero<-data_baseline$spc+abs(yminset)
#data_baseline_zero_hyperSpec<-new ("hyperSpec",
                                   #data = data.frame (Cells_bgsub[,1:ncol_meta]),
                                   #spc = data_baseline_zero,
                                   #wavelength=as.numeric(colnames(all_spc)[2:length(colnames(all_spc))]))
#plot (data_baseline_zero_hyperSpec, col = 1 : 2)
#write.csv(data_baseline_zero_hyperSpec,paste(output,"Cells_bg_baseline_zero.csv",sep=""),quote = F,row.names = F)

#################
##normalization##
#################
#data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowMeans (data_baseline_zero_hyperSpec)
data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowSums (data_baseline_zero_hyperSpec)
#normalization

data_baseline_normalization <- data_baseline_zero_scale_hyperSpec
#typeof(wave_num)
wave_num<-as.numeric(t(wave_num))
data_baseline_normalization <- Spc_at_wavenumber(data_baseline_normalization$spc,wave_num)
data_meta <- as.data.frame(data_baseline_zero_scale_hyperSpec@data[,-ncol(data_baseline_zero_scale_hyperSpec@data)])

data_meta <- select(data_meta, c(Name_group, "Group"))

data_baseline_normalization <- cbind(data_baseline_normalization, data_meta, RBCS_group=Meta_data$RBCS_group)
hyperspec_melt_summary <- Mean_spc(data_baseline_normalization,c(Name_group, "Group"))
hyperspec_melt_summary
RBCS_selectpeak_res <- RBCS_selectpeak(hyperspec_melt_summary, Pvalue = 0.0001)
write.table(RBCS_selectpeak_res,paste(outpath,"/result.txt",sep=''),sep = '\t',quote = F,col.names = T,row.names = T)

