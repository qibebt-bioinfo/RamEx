#################################################################
# Function:  Ramanome Pretreatment function
#Author: Rongzhe Chen, Yuehui He  
# Last update: 2021-11-01, Rongze Chen

#################################################################
# install necessary libraries

## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
sourcedir <- paste(sourcedir,"/databases/QC/Function")
source(paste(sourcedir, "util_ramex.R", sep = "/")) 
source(paste(sourcedir, "baseline.R", sep = "/"))
source(paste(sourcedir, "bg.R", sep = "/"))
source(paste(sourcedir, "CAST-R_function.R", sep = "/"))
source(paste(sourcedir, "CD-ratio.R", sep = "/"))
source(paste(sourcedir, "depend_packages.R", sep = "/"))
source(paste(sourcedir, "Draw.R", sep = "/"))
#source(paste(sourcedir, "Draw_spc_all.R", sep = "/"))
source(paste(sourcedir, "folder_reader.R", sep = "/"))
source(paste(sourcedir, "getspc.R", sep = "/"))
source(paste(sourcedir, "hclust_hyperspec.R", sep = "/"))
source(paste(sourcedir, "hyperspec_to_dataframe.R", sep = "/"))
source(paste(sourcedir, "initialization.R", sep = "/"))
source(paste(sourcedir, "install_packages_for_me.R", sep = "/"))
source(paste(sourcedir, "Mean_spc.R", sep = "/"))
source(paste(sourcedir, "meta_data.R", sep = "/"))
source(paste(sourcedir, "new_hyperSpec.R", sep = "/"))
source(paste(sourcedir, "Normalization.R", sep = "/"))
source(paste(sourcedir, "PCA.R", sep = "/"))
source(paste(sourcedir, "RamanEx_plot.R", sep = "/"))
source(paste(sourcedir, "Draw_spc_all.R", sep = "/"))
source(paste(sourcedir, "Real_time_box_plot.R", sep = "/"))
source(paste(sourcedir, "Real_time_calculation_system_Pvalue.R", sep = "/"))
source(paste(sourcedir, "Real_time_calculation_system.R", sep = "/"))
source(paste(sourcedir, "Real-time monitoring of minimum sample size.R", sep = "/"))
source(paste(sourcedir, "Renishaw_data_reader.R", sep = "/"))
source(paste(sourcedir, "Rtsne.R", sep = "/"))
source(paste(sourcedir, "Samplesize_plot.R", sep = "/"))
source(paste(sourcedir, "Sampling_depth2.R", sep = "/"))
source(paste(sourcedir, "Save_hyperspec.R", sep = "/"))
source(paste(sourcedir, "S-G_smooth.R", sep = "/"))
source(paste(sourcedir, "Single_cell_folder.R", sep = "/"))
source(paste(sourcedir, "Single-Spec_by group_with shadow.R", sep = "/"))
source(paste(sourcedir, "SNR_function.R", sep = "/"))
source(paste(sourcedir, "spc_melt.R", sep = "/"))
source(paste(sourcedir, "spc_quanlity_control.R", sep = "/"))
source(paste(sourcedir, "spc_waveumber.R", sep = "/"))
source(paste(sourcedir, "tSNE_plot.R", sep = "/"))

args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--raman_data_folder"),type="character", help="Input_the_raman_data_floder"),
        make_option(c("-f", "--data_format"), type="character", default="Horiba", help="Input_the_data_format[default %default]"),
	make_option(c("-p", "--perplexity"), type="character", default="5", help="Input_the_perplexity_of_the_tSNE[default %default]"),
	make_option(c("-o", "--out_dir"), type="character", default='Baseline_result', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$raman_data_folder)) stop('Please input a test raman file')

ramanpath <- opts$raman_data_folder 
outpath <- opts$out_dir
data_format <- opts$data_format
perplexity <- as.numeric(opts$perplexity)

#outputpath creation
dir.create(outpath)


#读取数据并根据命名规则生成metadata文件 
Name_group <- c("group_A","group_B","group_C")  

###数据标准化
if(data_format != "Horiba"){
	data_group <- renishaw(ramanpath,Name_group,save = T,outpath) # renishaw的数据
	fwrite(data_group,paste(outpath,"wavenumber_group.txt",sep = "/"),row.names = F,col.names = T,quote = F,sep = "\t")
	all_spc <- folder_reader(paste(outpath,"all_single_cell",sep = "/")) } else{
	all_spc <- folder_reader(ramanpath)  
}
###平滑
#all_data <- S_G_smooth(all_spc)
data_hyperSpec <- new_hyperSpec(all_spc,Name_group)
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
#data_baseline_normalization <- Normalization(data_hyperSpec,from = 2500,to = 3000,nor_cal = "max")
data_baseline_normalization_plot <- Spec_all(data_baseline_normalization)
#fwrite(data_baseline_normalization_plot, paste(outpath,"All_data.txt",sep = "/"), sep = "\t")
ggsave(filename=paste(outpath,"/","data_baseline_nor_all.png", sep=""),plot=data_baseline_normalization_plot, limitsize=T,width = 16,height = 4)
data_baseline_normalization_all <- cbind(Point=data_baseline_normalization$filename,data.frame(data_baseline_normalization$spc))
#fwrite(data_baseline_normalization_plot,paste(outpath,"All_data.txt",sep = "/"),row.names = F,col.names = T,quote = F,sep = "\t")
data_baseline_meta <- data.frame(data_baseline_normalization@data[,-ncol(data_baseline_normalization@data)])
write.table(data_baseline_meta, file=paste(outpath,"All_Meta.txt",sep = "/"), append = FALSE, row.name= FALSE, quote = FALSE, sep = "\t")
write.table(data_baseline_normalization_all, file=paste(outpath,"All_data.txt",sep = "/"), append = FALSE, row.name= FALSE, quote = FALSE, sep = "\t")
####t-SNE
tSNE_data <- tSNE(data_baseline_normalization$spc,data_baseline_normalization$Group,outpath = outpath,perplexity= perplexity)
tSNE_data_plot <- tsne_plot(tSNE_data)
ggsave(filename=paste(outpath,"/","tsne.png", sep=""),plot=tSNE_data_plot, limitsize=T,width = 6,height = 4)

###PCA
data_pca <- pca(data_baseline_normalization$spc,data_baseline_normalization$Group)
data_pca_plot <- pca_plot(data_pca, data_baseline_normalization$Group)
ggsave(filename=paste(outpath,"/","data_pca.png", sep=""),plot=data_pca_plot, limitsize=T,width = 6,height = 4)

###PLS-DA
plsda.res <- PLSDA(data_baseline_normalization$spc,data_baseline_normalization$Group)
plsda.res_plot <- plsda_plot(plsda.res,name = "plsda",data_baseline_normalization$Group)
ggsave(filename=paste(outpath,"/","data_plsda.png", sep=""),plot=plsda.res_plot, limitsize=T,width = 6,height = 4)

###
