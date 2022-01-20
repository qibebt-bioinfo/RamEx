
#################################################################
# Function:  Ramanome IRCA Generate Local network by select bands function
# Author: Yuehui He, Shi Huang
# Last update: 2021-09-24, Yuehui He, Shi Huang
#################################################################
# install necessary libraries
library(factoextra)
## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util_clean.R',sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--clean_data"),type="character", help="Input_the_raman_clean_data"),
	make_option(c("-m", "--meta_data"),type="character", help="Input_the_meta_data"),
	make_option(c("-o", "--out_dir"), type="character", default='local_network', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$clean_data)) stop('Please input the clean raman data')
if(is.null(opts$meta_data)) stop('Please input the meta data')
matrixfile<-opts$clean_data
metadatafile<-opts$meta_data
plot_local_IRCN<-TRUE
LocalBandsAnnfile<- sprintf('%s/databases/Local_bands_annotation.txt',sourcedir) 
if(!is.null(LocalBandsAnnfile)){ local_bands_ann<-read.table(LocalBandsAnnfile,header=T, sep="\t") }

outpath <- opts$out_dir#"outputpath"
category<-c("group_A", "group_B", "group_C")



#outputpath creation
dir.create(outpath)
outpath_png <- paste(outpath,"/Local-IRCA_png/",sep="")
dir.create(outpath_png)
options(warn=-1)
#-------------------------------
# Metadata input
#-------------------------------
allmetadata<-read.table(metadatafile,header=T,sep="\t",row.names=1); 
metadata<-data.frame(allmetadata[order(rownames(allmetadata)), ]) 
colnames(metadata)<-colnames(allmetadata)
all_group<-colnames(metadata)
all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]; try(metadata_f<-metadata[, all_group_f])
all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]; try(metadata_n<-metadata[, all_group_n])

if(length(category)>1){
  metadata$group_c <- as.factor(apply( metadata[ , category ], 1 , paste, collapse = "__" ))
  ori_category<-category
  category<-paste(ori_category, collapse="__")
  names(metadata)[length(metadata)]<-category
  all_group_f[length(all_group_f)+1]<-category
}else{
  ori_category<-category
}
Group<-metadata[, category]
#-------------------------------
# cleandata input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t"); 
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-formatC(raw_wn, digits=0, format="f")
cor_mats_rpvalue<-by(as.matrix(mat),Group, rcorr_df)
v_cor_mats_rpvalue<-lapply(cor_mats_rpvalue, function(x) vectorize_dm_rcorr(x, group=NULL, duplicate=FALSE))
#-------------------------------
# Local IRCNs by chordDiagram
# Generate local networks by selected bands
#-------------------------------

Peak_sub<-as.character(local_bands_ann$Wave_num)
bands<-Peak_sub
Peaks_fixed_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
names(Peaks_fixed_cor_mats_rpvalue)<-names(cor_mats_rpvalue)
Peak_sub<-local_bands_ann$Wave_num
Peak_sub
#-------------------------------
# Generate local networks by selected bands
#-------------------------------
bands<-as.character(Peak_sub)
bands_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
local_bands_ann_new<-local_bands_ann
local_bands_ann_new$Wave_num<-as.character(local_bands_ann_new$Wave_num)
local_bands_ann_new$Wave_num<-as.character(round(as.numeric(gsub("B","",as.character(local_bands_ann_new$Wave_num))),0))

#---------------------------------------------------
# png plot_local_chorddiagram
#---------------------------------------------------
# 1.# Pos_Edge=F, Neg_Edge=T, Threshold=0.6,
#---------------------------------------------------
device="png"
outpath_local_neg6_png_rpvalue <-paste(outpath_png,"/local_neg6_png_rpvalue/",sep = "")
dir.create(outpath_local_neg6_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
	Cairo(file = paste(outpath_local_neg6_png_rpvalue, "Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
	Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0.6, local_bands_ann_new)
	dev.off()  
}
#-----------------------------------------------------------------------------------------------
#---------------------------------------------------
# 2.# Pos_Edge=T, Neg_Edge=F, Threshold=0.6,
#---------------------------------------------------
device="png"
outpath_local_pos6_png_rpvalue <- paste(outpath_png,"/local_pos6_png_rpvalue/",sep = "")
dir.create(outpath_local_pos6_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
	Cairo(file = paste(outpath_local_pos6_png_rpvalue, "Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0.6, local_bands_ann_new)
    dev.off()  
}
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 3.# Pos_Edge=T, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
device="png"
outpath_local_neg6pos6_png_rpvalue <- paste(outpath_png,"/local_neg6pos6_png_rpvalue/",sep = "")
dir.create(outpath_local_neg6pos6_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
  Cairo(file = paste(outpath_local_neg6pos6_png_rpvalue, "Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
  Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0.6, local_bands_ann_new)
  dev.off()  
}
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 4.# Pos_Edge=F, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
device="png"
outpath_local_neg0_png_rpvalue <-  paste(outpath_png,"local_neg0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0_png_rpvalue/"
dir.create(outpath_local_neg0_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
	Cairo(file = paste(outpath_local_neg0_png_rpvalue, "Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0, local_bands_ann_new)
    dev.off()  
}
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 5.# Pos_Edge=T, Neg_Edge=F, Threshold=0,
  #---------------------------------------------------
device="png"
outpath_local_pos0_png_rpvalue <- paste(outpath_png,"local_pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_pos0_png_rpvalue/"
dir.create(outpath_local_pos0_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
  Cairo(file = paste(outpath_local_pos0_png_rpvalue, "Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
  Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0, local_bands_ann_new)
  dev.off()  
}
#-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 6.# Pos_Edge=T, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
device="png"
outpath_local_neg0pos0_png_rpvalue <- paste(outpath_png,"/local_neg0pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0pos0_png_rpvalue/"
dir.create(outpath_local_neg0pos0_png_rpvalue)
for (group_name in names(bands_cor_mats_rpvalue)) {
  Cairo(file = paste(outpath_local_neg0pos0_png_rpvalue, "/Local_chordDiagram_", group_name, ".",device , sep=""), unit="in", dpi=300, width=30, height=30, type=device, bg="white")
   Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0, local_bands_ann_new)
    dev.off()  
}
  #----------------------------------------------------------------------------------------------
