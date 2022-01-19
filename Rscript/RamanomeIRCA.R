
#################################################################
# Function:  Ramanome IRCA analysis  
# Last update: 2021-08-18, He Yuehui; Shi Huang
#################################################################

# install necessary libraries


p <- c("optparse","RColorBrewer","igraph", "circlize", "permute", "ggplot2", "reshape2", "proxy", "pheatmap", "Cairo", "grid", "plyr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
	#install.packages(p, dep=TRUE)<F12>
    suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- Sys.getenv("RamEx")
source(sprintf('%s/Rscript/util.R',sourcedir))


args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
	make_option(c("-i", "--raman_data"),type="character", help="Input_the_raman_data"),
	make_option(c("-m", "--meta_data"), type="character", help="Input_raman_meta_data"),
	make_option(c("-o", "--out_dir"), type="character", default='Ramandata', help="outpath file name[default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)
if(is.null(opts$raman_data)) stop('Please input a test raman file')
if(is.null(opts$meta_data)) stop('Please input a background data file')

matrixfile <- opts$raman_data # "approxfun_Alldata.txt"
metadatafile <- opts$meta_data #"approxfun_Alldata_Map.txt"

#---------------------
# parameter input
#---------------------

Pos_Edge<-FALSE
Neg_Edge<-TRUE
GlobalBandsAnnfile<- sprintf('%s/databases/Global_bands_annotation.txt',"/home/gene/jinggc/RamEX") #GlobalBandsAnnfile<-NULL
LocalBandsAnnfile<- sprintf('%s/databases/Local_bands_annotation.txt',"/home/gene/jinggc/RamEX") 
category<-c("Celltype","Timepoint")#c("Species","Strain","Condition","Timepoint","User","Date")
cor_cutoff<-0.6

plot_local_IRCN<-TRUE
plot_global_IRCN<-TRUE
#device="pdf"

#---------------------
# output pathway
#---------------------
outpath <- "./Timepoint_neg0.6/"
dir.create(outpath)
outpath_png <- paste(outpath,"Global-Local-IRCA_png/",sep="")
dir.create(outpath_png)

#-------------------------------
# Spec range trimming
#-------------------------------
#SpecRange_remove<-1800:2600

options(warn=-1)
#-------------------------------
# Add prefix to outpath
#-------------------------------
#prefix<-paste(paste(ifelse(Pos_Edge, "Pos", "NoPos"),ifelse(Neg_Edge, "Neg", "NoNeg"),sep="-"), cor_cutoff, "", sep="_")
#outpath<-paste(outpath, prefix, sep="")
#-------------------------------
# Spectral data input
#-------------------------------
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t"); 
#remove 3d_s1-s2outliers
#mat<-mat[-c(603,610,623,637,631),]
#mat<-mat[order(rownames(mat)),]

#-------------------------------
# scale by dividing sum area
#-------------------------------
scale_sum<-function(data){
  sum<-sum(data)
  nor<-data/sum
}
mat<-apply(mat, 1, scale_sum)
mat<-data.frame(t(mat))
#write.table(mat,paste(outpath,"Nor_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
write.csv(mat,paste(outpath,"Nor_",matrixfile,".csv",sep=""),row.names=T,quote=F)
#mat_show<-mat[1:5,1:5]
#-------------------------------
# Clean data (Negative, NA, Inf, or -Inf)
#-------------------------------
mat<-CleanData(mat)
cat("The number of Raman shifts: ", ncol(mat) , "\n")
#-------------------------------
# Wave number trimming
#-------------------------------
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-formatC(raw_wn, digits=0, format="f")
if(length(unique(wn0)) == length(raw_wn)){
  colnames(mat)<-paste("B", wn0,sep="")
}else{
  wn1<-formatC(raw_wn, digits=1, format="f")
  colnames(mat)<-paste("B", wn1, sep="")
}
#-------------------------------
# Spec range trimming
#-------------------------------
#SpecRange_remove<-1800:2600
#SpecCols_remove<-paste("B", SpecRange_remove, sep="")
#mat<-mat[, which(!colnames(mat)%in% SpecCols_remove)]
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-round(raw_wn, 0)
cat("The number of Raman shifts: ", ncol(mat) , "\n")
#-------------------------------
# Metadata input
#-------------------------------
allmetadata<-read.table(metadatafile,header=T,sep="\t",row.names=1); 
#remove 3d_s1-s2outliers
#allmetadata<-allmetadata[-c(603,610,623,637,631),]
#metadata<-data.frame(allmetadata[order(rownames(allmetadata)), which(sapply(allmetadata,var)!=0)]) 
#colnames(metadata)<-colnames(allmetadata)[which(sapply(allmetadata,var)!=0)]
metadata<-data.frame(allmetadata[order(rownames(allmetadata)), ]) 
colnames(metadata)<-colnames(allmetadata)
all_group<-colnames(metadata)
all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]; try(metadata_f<-metadata[, all_group_f])
all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]; try(metadata_n<-metadata[, all_group_n])

#-------------------------------
# Bands annotation file input
#-------------------------------
if(!is.null(GlobalBandsAnnfile)){ global_bands_ann<-read.table(GlobalBandsAnnfile,header=T, sep="\t") 
Wave_num <- data.frame(Wave_num = wn0)
global_bands_ann <- merge(Wave_num, global_bands_ann, by="Wave_num", all.x=TRUE, sort = TRUE)
global_bands_ann[is.na(global_bands_ann)] <- "Unknown"
global_bands_ann$Wave_num<-as.numeric(global_bands_ann$Wave_num)
global_bands_ann<-global_bands_ann[with(global_bands_ann, order(Wave_num)),]

}else{
  stop("Please provide global Bands annotation file!")
}
if(!is.null(LocalBandsAnnfile)){ local_bands_ann<-read.table(LocalBandsAnnfile,header=T, sep="\t") }

#-------------------------------
# If multiple categories specified
#-------------------------------
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
# Generate correlation matrix by one category in metadata
#-------------------------------
cor_mats<-by(mat, Group, cor)
if(!is.null(GlobalBandsAnnfile)){
                                 v_cor_mats<-lapply(cor_mats, function(x) vectorize_dm(x, group=global_bands_ann$Group, duplicate=FALSE) )
                                 }else{
                                 v_cor_mats<-lapply(cor_mats, function(x) vectorize_dm(x, group=NULL, duplicate=FALSE)) 
                                 }
#write.csv(v_cor_mats["CC4324__00h"],paste(outpath,"CC4324_00h.csv"))
v_cor_mats_trimmed<-lapply(v_cor_mats, function(x) Trim_v_corr_mat(x, Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge, Threshold=cor_cutoff))
#write.csv(v_cor_mats_trimmed[1],paste(outpath,"1_v_cor_mats_trimmed.csv"))


#-------------------------------
# Global network by igraph
#-------------------------------
g_array<-sapply(names(cor_mats), function(x) Plot_network_graph(cor_mats[[x]], Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge,Threshold=cor_cutoff, node_size=1, node_color="Cluster_membership", layout="layout_components", outdir=NULL)) # paste(outpath, x, "_global_network_igraph.pdf",sep="")
g<-g_array["g",]
#g_array<-sapply(names(cor_mats), function(x) Plot_network_graph(cor_mats[[x]], Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge,Threshold=cor_cutoff, node_size=1, node_color="Cluster_membership", layout="layout_components", outdir=paste(outpath, x, "_global_network_igraph.pdf",sep=""))) # paste(outpath, x, "_global_network_igraph.pdf",sep="")
#g<-g_array["g",]

#-------------------------------
# Global network statistics by igraph
#-------------------------------
out_array<-t(g_array[1:16,])
net_stats<-data.frame(matrix(unlist(out_array), nrow=nrow(out_array), byrow=F))
dimnames(net_stats)<-dimnames(out_array)
sink(paste(outpath,"global_network_stats_by_", category,".xls",sep="")); cat("\t"); write.table(net_stats,sep="\t",quote=FALSE); sink()
#-------------------------------
# Visualization: Global network statistics by ggplot2
#-------------------------------
net_stats_df <-melt(data.frame( Group=rownames(net_stats), net_stats))
p<-ggplot(net_stats_df, aes(x=Group, y=value))+ geom_bar(stat = "identity") + 
   xlab( category ) + ylab("Network-level topological feature") + 
   coord_flip()+
   facet_grid(. ~ variable, scales="free") 
ggsave(filename=paste(outpath, "global_network_stats_barplot_by_", category,".ggplot.png", sep=""), plot=p, limitsize=FALSE, width=24, height=1 + nlevels(Group)*0.2, device='png') 

#-------------------------------
# Cluster membership distribution
#-------------------------------
Cluster_membership<-sapply(g, function(x) clusters(x, mode = "strong")$membership)
Cluster_membership_df<-data.frame(Cluster_membership, global_bands_ann)
#Cluster_membership_df<-data.frame(Cluster_membership)
sink(paste(outpath,"global_network_Cluster_membership_by_", category,".xls",sep="")); cat("\t"); write.table(Cluster_membership_df,sep="\t",quote=FALSE); sink()

Cluster_walktrap_membership<-sapply(g, function(x) cluster_walktrap(x)$membership)
Cluster_walktrap_membership_df<-data.frame(Cluster_walktrap_membership, global_bands_ann)
sink(paste(outpath,"global_network_Cluster_walktrap_membership_by_", category,".xls",sep="")); cat("\t"); write.table(Cluster_walktrap_membership_df,sep="\t",quote=FALSE); sink()

#-------------------------------
# Degree distribution
#-------------------------------
Degree<-sapply(g, function(x) igraph::degree(x)); 
Degree_c<-data.frame(Wave_num=as.numeric(gsub("[A-Z]", "", rownames(Degree))), Degree)
Degree_df<-merge(global_bands_ann, Degree_c, by=1); 
rownames(Degree_df)<-paste("B", Degree_df$Wave_num, sep="")
sink(paste(outpath,"global_network_Degree_by_", category,".xls",sep=""));
cat("\t"); write.table(Degree_df,sep="\t",quote=FALSE); sink()
#-------------------------------
    d<-Degree_df[, -(1:3)]
    ann = data.frame(Degree_df[, 2])
    rownames(d)<-rownames(ann)<-paste(Degree_df[, 1], Degree_df[, 3], sep="__")
    ann_colors = list(Group = ann[, 1])
#    h<-pheatmap(t(d),annotation =ann, cluster_row=F,cluster_col=F, fontsize_row=4, fontsize_col=2, cellheight=4, cellwidth=2, fontsize=4, filename=paste(outpath, "Degree.heatmap.pdf",sep=""))


#--------------------------------------------------------
# Generate correlation matrix by one category in metadata
# calculate both correlation and p-value
# vectorize both correlation and p-value
#--------------------------------------------------------
require(devtools)
#install_version("data.table", version = "1.11.8", repos = "http://cran.us.r-project.org")
#install.packages("/mnt/data6/heyh/Network_20200514/Hmisc_4.4-0.tar.gz")
#source("/mnt/data6/heyh/Network_20200514/Function_v_cor_mat_rpvalue_util.R")
cor_mats_rpvalue<-by(as.matrix(mat),Group, rcorr_df)
v_cor_mats_rpvalue<-lapply(cor_mats_rpvalue, function(x) vectorize_dm_rcorr(x, group=NULL, duplicate=FALSE)) 
#save(cor_mats_rpvalue,v_cor_mats_rpvalue,file = paste(outpath, "v_cor_mats_rpvalue.RData",sep="")) 

#---------------------------------------------------------------------------------
#cor_msts transfer to cor_mats_rpvalue for Plot_global_chordDiagram_rpvalue
#---------------------------------------------------------------------------------
cor_mats_rpvalue_sig005<-NULL
for (list_name in names(cor_mats_rpvalue)) {
  #list_name<-"000h"
  cor_mats_r_tem<-cor_mats_rpvalue[[list_name]]$r
  cor_mats_pvalue_tem<-cor_mats_rpvalue[[list_name]]$P
  cor_mats_r_tem[which(cor_mats_pvalue_tem>0.05)]<-0 # p.value>0.05
  #cor_mats_r_tem[which(cor_mats_r_tem>(-0.6))]<-0 # p.value>0.05
  cor_mats_rpvalue_tem<-list(cor_mats_r_tem)
  names(cor_mats_rpvalue_tem)<-list_name
  cor_mats_rpvalue_sig005<-c(cor_mats_rpvalue_sig005,cor_mats_rpvalue_tem)
}
#cor_mats_rpvalue_sig005<-cor_mats_rpvalue_new

#-------------------------------------------------------------------------------------------------------
# calculate connectedness P<0.05 (i.e.MeansNegCorr / MeansPosCorr / SumsNegCorr / SumsPosCorr/ Degree_Neg / Degree_Pos)
# cor_mats_rpvalue_new (P<0.05) cor_cutoff<-0
#--------------------------------------------------------------------------------------------------------
outpath_global_stat<-paste(outpath,"global_stat/",sep="")
dir.create(outpath_global_stat)
#save(cor_mats_rpvalue_sig005,paste(outpath_global_stat,"cor_mats_rpvalue_sig005.RData",sep=""))

cor_cutoff_new<-0
SumsNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
SumsPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
MeansNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
MeansPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
CountsNegCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
CountsPosCorr_cutoff0<-data.frame(A=rep("A",ncol(mat)))
for (group in names(cor_mats_rpvalue_sig005)) {
  #group="000h"
  a1<-cor_mats_rpvalue_sig005[[group]]
  a2<-a1
  a2[a2<cor_cutoff_new]<-0 #postive corr
  a3<-a1
  a3[a3>(-cor_cutoff_new)]<-0 #negative corr
  s2<-colSums(a2)#postive corr
  s3<-colSums(a3)#negative corr
  SumsNegCorr_cutoff0<-cbind(SumsNegCorr_cutoff0,data.frame(s3))
  colnames(SumsNegCorr_cutoff0)[ncol(SumsNegCorr_cutoff0)]<-group
  SumsPosCorr_cutoff0<-cbind(SumsPosCorr_cutoff0,data.frame(s2))
  colnames(SumsPosCorr_cutoff0)[ncol(SumsPosCorr_cutoff0)]<-group
  
  s4<-colMeans(a2)#postive corr
  s5<-colMeans(a3)#negative corr
  MeansNegCorr_cutoff0<-cbind(MeansNegCorr_cutoff0,data.frame(s5))
  colnames(MeansNegCorr_cutoff0)[ncol(MeansNegCorr_cutoff0)]<-group
  MeansPosCorr_cutoff0<-cbind(MeansPosCorr_cutoff0,data.frame(s4))
  colnames(MeansPosCorr_cutoff0)[ncol(MeansPosCorr_cutoff0)]<-group
  
  a2[a2<=0]<-0 
  a2[a2>0]<-1 #count numbers of postive corr 
  a3[a3>=0]<-0 
  a3[a3<0]<-1 #count numbers of negative corr
  count2<-colSums(a2)
  count3<-colSums(a3)
  CountsNegCorr_cutoff0<-cbind(CountsNegCorr_cutoff0,data.frame(count3))
  colnames(CountsNegCorr_cutoff0)[ncol(CountsNegCorr_cutoff0)]<-group
  CountsPosCorr_cutoff0<-cbind(CountsPosCorr_cutoff0,data.frame(count2))
  colnames(CountsPosCorr_cutoff0)[ncol(CountsPosCorr_cutoff0)]<-group
}
SumsNegCorr_cutoff0<-SumsNegCorr_cutoff0[-1]
SumsPosCorr_cutoff0<-SumsPosCorr_cutoff0[-1]
MeansNegCorr_cutoff0<-MeansNegCorr_cutoff0[-1]
MeansPosCorr_cutoff0<-MeansPosCorr_cutoff0[-1]
CountsNegCorr_cutoff0<-CountsNegCorr_cutoff0[-1]
CountsPosCorr_cutoff0<-CountsPosCorr_cutoff0[-1]
#write.csv(SumsNegCorr,paste(outpath,"rpvalue_SumsNegCorr_by_",category,".csv",sep=""))
#write.csv(SumsPosCorr,paste(outpath,"rpvalue_SumsPosCorr_by_",category,".csv",sep=""))
#write.csv(MeansNegCorr,paste(outpath,"rpvalue_MeansNegCorr_by_",category,".csv",sep=""))
#write.csv(MeansPosCorr,paste(outpath,"rpvalue_MeansPosCorr_by_",category,".csv",sep=""))
#write.csv(CountsNegCorr,paste(outpath,"rpvalue_Degree_Neg_cutoff_",cor_cutoff,"_CountsNegCorr_by_",category,".csv",sep=""))
#write.csv(CountsPosCorr,paste(outpath,"rpvalue_Degree_Pos_cutoff_",cor_cutoff,"_CountsPosCorr_by_",category,".csv",sep=""))
#write.csv(CountsNegCorr_cutoff0,paste(outpath,"rpvalue_Degree_Neg_cutoff_",0,"_CountsNegCorr_by_",category,".csv",sep=""))
#write.csv(CountsPosCorr_cutoff0,paste(outpath,"rpvalue_Degree_Pos_cutoff_",0,"_CountsPosCorr_by_",category,".csv",sep=""))
write.csv(SumsNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_SumsNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(SumsPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_SumsPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_MeansNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_MeansPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(CountsNegCorr_cutoff0,paste(outpath_global_stat,"rpvalue_Degree_Neg_cutoff_",cor_cutoff_new,"_CountsNegCorr.csv",sep=""))
write.csv(CountsPosCorr_cutoff0,paste(outpath_global_stat,"rpvalue_Degree_Pos_cutoff_",cor_cutoff_new,"_CountsPosCorr.csv",sep=""))
#--------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# calculate connectedness P<0.05 (i.e.MeansNegCorr / MeansPosCorr / SumsNegCorr / SumsPosCorr/ Degree_Neg / Degree_Pos)
# cor_mats_rpvalue_new (P<0.05)
#--------------------------------------------------------------------------------------------------------
cor_cutoff_new<-0.6
SumsNegCorr<-data.frame(A=rep("A",ncol(mat)))
SumsPosCorr<-data.frame(A=rep("A",ncol(mat)))
MeansNegCorr<-data.frame(A=rep("A",ncol(mat)))
MeansPosCorr<-data.frame(A=rep("A",ncol(mat)))
CountsNegCorr<-data.frame(A=rep("A",ncol(mat)))
CountsPosCorr<-data.frame(A=rep("A",ncol(mat)))
for (group in names(cor_mats_rpvalue_sig005)) {
  #group="000h"
  a1<-cor_mats_rpvalue_sig005[[group]]
  a2<-a1
  a2[a2<cor_cutoff_new]<-0 #postive corr
  a3<-a1
  a3[a3>(-cor_cutoff_new)]<-0 #negative corr
  s2<-colSums(a2)#postive corr
  s3<-colSums(a3)#negative corr
  SumsNegCorr<-cbind(SumsNegCorr,data.frame(s3))
  colnames(SumsNegCorr)[ncol(SumsNegCorr)]<-group
  SumsPosCorr<-cbind(SumsPosCorr,data.frame(s2))
  colnames(SumsPosCorr)[ncol(SumsPosCorr)]<-group
  
  s4<-colMeans(a2)#postive corr
  s5<-colMeans(a3)#negative corr
  MeansNegCorr<-cbind(MeansNegCorr,data.frame(s5))
  colnames(MeansNegCorr)[ncol(MeansNegCorr)]<-group
  MeansPosCorr<-cbind(MeansPosCorr,data.frame(s4))
  colnames(MeansPosCorr)[ncol(MeansPosCorr)]<-group
  
  a2[a2<cor_cutoff_new]<-0 
  a2[a2>=cor_cutoff_new]<-1 #count numbers of postive corr 
  a3[a3>(-cor_cutoff_new)]<-0 
  a3[a3<=(-cor_cutoff_new)]<-1 #count numbers of negative corr
  count2<-colSums(a2)
  count3<-colSums(a3)
  CountsNegCorr<-cbind(CountsNegCorr,data.frame(count3))
  colnames(CountsNegCorr)[ncol(CountsNegCorr)]<-group
  CountsPosCorr<-cbind(CountsPosCorr,data.frame(count2))
  colnames(CountsPosCorr)[ncol(CountsPosCorr)]<-group
}
SumsNegCorr<-SumsNegCorr[-1]
SumsPosCorr<-SumsPosCorr[-1]
MeansNegCorr<-MeansNegCorr[-1]
MeansPosCorr<-MeansPosCorr[-1]
CountsNegCorr<-CountsNegCorr[-1]
CountsPosCorr<-CountsPosCorr[-1]
write.csv(SumsNegCorr,paste(outpath_global_stat,"rpvalue_SumsNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(SumsPosCorr,paste(outpath_global_stat,"rpvalue_SumsPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansNegCorr,paste(outpath_global_stat,"rpvalue_MeansNegCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(MeansPosCorr,paste(outpath_global_stat,"rpvalue_MeansPosCorr_cutoff_",cor_cutoff_new,".csv",sep=""))
write.csv(CountsNegCorr,paste(outpath_global_stat,"rpvalue_Degree_Neg_cutoff_",cor_cutoff_new,"_CountsNegCorr.csv",sep=""))
write.csv(CountsPosCorr,paste(outpath_global_stat,"rpvalue_Degree_Pos_cutoff_",cor_cutoff_new,"_CountsPosCorr.csv",sep=""))
#--------------------------------------------------------------------------------------------------------


#load(paste("/mnt/data6/heyh/IRCN_CC124_20200624/Timepoint_neg0.6/","mydata_by_Celltype__Timepoint.RData",sep=""))
#------------------------------------------------------------------
# Global IRCNs
# E:\RWAS\Plot_Fig8BC_chordDiagram_20200623\Plot_Fig8BC_global_chordDiagram_20200624.R
#------------------------------------------------------------------
if(plot_global_IRCN){
  #---------------------------------------------------
  # 1.1 # Pos_Edge=F, Neg_Edge=T, Threshold=0.6, Degree
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_global_neg6_png_rpvalue <-paste(outpath_png,"global_neg6_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_neg6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N-__000h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsNegCorr))),
                         y=CountsNegCorr[,group_name]/max(CountsNegCorr))
    #CairoPNG(file = paste(outpath_global_neg6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png", sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_neg6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=FALSE, 
                                     Neg_Edge=TRUE, 
                                     Threshold=0.6)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 2.1 # Pos_Edge=T, Neg_Edge=F, Threshold=0.6, degree
  #---------------------------------------------------
  #source("/mnt/data6/heyh/Network_20200426/Plot_local_chordDiagram_new_20200426.r")
  #load(file = paste(outpath,"NoPos-Neg_0.8_mydata_by_Timepoint.RData",sep=""))
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_global_pos6_png_rpvalue <-paste(outpath_png,"global_pos6_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_pos6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N-__000h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
	#group_name<-"CC124_072h"
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsPosCorr))),
                         y=CountsPosCorr[,group_name]/max(CountsPosCorr))
    #CairoPNG(file = paste(outpath_global_pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png" , sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=TRUE, 
                                     Neg_Edge=FALSE, 
                                     Threshold=0.6)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 3.# Pos_Edge=T, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  #source("/mnt/data6/heyh/Network_20200426/Plot_local_chordDiagram_new_20200426.r")
  #load(file = paste(outpath,"NoPos-Neg_0.8_mydata_by_Timepoint.RData",sep=""))
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_global_neg6pos6_png_rpvalue <- paste(outpath_png,"global_neg6pos6_png_rpvalue/",sep="")
  dir.create(outpath_global_neg6pos6_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    #  a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(CountsPosCorr))),
    #                    y=CountsPosCorr[,group_name]/max(CountsPosCorr))
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    #CairoPNG(file = paste(outpath_global_neg6pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png" , sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_neg6pos6_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=NULL,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=TRUE, 
                                     Neg_Edge=TRUE,
                                     Threshold=0.6)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 4.1 # Pos_Edge=F, Neg_Edge=T, Threshold=0,degree
  #---------------------------------------------------
  #source("/mnt/data6/heyh/Network_20200426/Plot_local_chordDiagram_new_20200426.r")
  #load(file = paste(outpath,"NoPos-Neg_0.8_mydata_by_Timepoint.RData",sep=""))
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_global_neg0_png_rpvalue <- paste(outpath_png,"global_neg0_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_neg0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsNegCorr_cutoff0))),
                         y=CountsNegCorr_cutoff0[,group_name]/max(CountsNegCorr_cutoff0))
    #CairoPNG(file = paste(outpath_global_neg0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png" , sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_neg0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=F, 
                                     Neg_Edge=T,
                                     Threshold=0)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  
  #---------------------------------------------------
  # 5.1 # Pos_Edge=T, Neg_Edge=F, Threshold=0,degree(CountsPosCorr_cutoff0)
  #---------------------------------------------------
  #source("/mnt/data6/heyh/Network_20200426/Plot_local_chordDiagram_new_20200426.r")
  #load(file = paste(outpath,"NoPos-Neg_0.8_mydata_by_Timepoint.RData",sep=""))
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_global_pos0_png_rpvalue <- paste(outpath_png,"global_pos0_png_rpvalue_degree/",sep="")
  dir.create(outpath_global_pos0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
                         x=as.numeric(gsub("[A-Z]", "", row.names(CountsPosCorr_cutoff0))),
                         y=CountsPosCorr_cutoff0[,group_name]/max(CountsPosCorr_cutoff0))
    #CairoPNG(file = paste(outpath_global_pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png", sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=a_degree,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=NULL,
                                     main_title=group_name,
                                     Pos_Edge=T, 
                                     Neg_Edge=F, 
                                     Threshold=0)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 6.# Pos_Edge=T, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  #source("/mnt/data6/heyh/Network_20200426/Plot_local_chordDiagram_new_20200426.r")
  #load(file = paste(outpath,"NoPos-Neg_0.8_mydata_by_Timepoint.RData",sep=""))
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_global_neg0pos0_png_rpvalue <- paste(outpath_png,"global_neg0pos0_png_rpvalue/",sep="")
  dir.create(outpath_global_neg0pos0_png_rpvalue)
  for (group_name in names(cor_mats_rpvalue_sig005)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    #  a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(SumsPosCorr_cutoff0))),
    #                     y=SumsPosCorr_cutoff0[,group_name]/max(SumsPosCorr_cutoff0))
    Corr_mat<-cor_mats_rpvalue_sig005[[group_name]]
    #CairoPNG(file = paste(outpath_global_neg0pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".png", sep=""), width=9, height=9, bg="white")
	Cairo(file = paste(outpath_global_neg0pos0_png_rpvalue, "Gobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=9, height=9, type=device, bg="white")
    #pdf(paste(outpath, "Global_chordDiagram_", group_name, ".pdf", sep=""), width=10, height=10)
    Plot_global_chordDiagram_rpvalue(Corr_mat=Corr_mat,
                                     a_degree=NULL,
                                     a_SumCorr=NULL,
                                     a_MeanCorr=a_MeanCorr,
                                     main_title=group_name,
                                     Pos_Edge=T,
                                     Neg_Edge=T,
                                     Threshold=0)
  }
  #-----------------------------------------------------------------------------------------------
}
#---------------------------------------------------------------------------------------------------------------------

# #--------------------------------------
# # Global network by pheatmap (plus HI, Degree, meanSCRS)
# #--------------------------------------
# #IRCN_neg6_anno
# #Degree_SCRS<-read.table(paste("E:/RWAS/IRCN_CC124_20200704/Timepoint_neg0.6/",
# #                              "global_network_Degree_by_Celltype__Timepoint.xls",sep = ""),
# #                        row.names = 1,header = T,
# #                        sep="\t")
# library(ComplexHeatmap)
# #Degree_SCRS<-Degree_df
# #Degree_SCRS<-t(Degree_SCRS[,-c(1,2,3)])
# Degree_SCRS<-t(Degree_df[,names(cor_mats_rpvalue)])
# Degree_SCRS_show<-Degree_SCRS[,1:5]
# outpath_heatmap_IRCN_neg6_anno<-paste(outpath,"Heatmap_IRCN_neg6_anno/",sep = "")
# dir.create(outpath_heatmap_IRCN_neg6_anno)
# for (name_mat in names(cor_mats_rpvalue)) {
#   #name_mat<-"CC124__168h"
#   #name_mat<-"CC124__000h"
#   Mean_SCRS<-apply(mat[which(Group==name_mat),],2,mean)
#   SD_SCRS<-apply(mat[which(Group==name_mat),],2,sd)
#   HI_SCRS<-SD_SCRS/Mean_SCRS
#   #ha_Mean_SCRS<-HeatmapAnnotation(foo = anno_lines(as.numeric(Mean_SCRS[1:50])), height = unit(2, "cm"))
#   #ha_HI_SCRS<-HeatmapAnnotation(foo = anno_lines(as.numeric(HI_SCRS[1:50])),height = unit(2, "cm"))
#   #ha_Degree_SCRS<-HeatmapAnnotation(foo = anno_lines(as.numeric(Degree_SCRS["CC124__168h",1:50])),height = unit(2, "cm"))
#   
#   mat_cor<-cor_mats_rpvalue[[name_mat]]
#   mat_cor_r<-data.matrix(mat_cor$r)
#   mat_cor_p<-data.matrix(mat_cor$P)
#   mat_cor_r[mat_cor_p>(0.05)]<-NA
#   mat_cor_r[upper.tri(mat_cor_r)]<-NA
#   mat_cor_r[mat_cor_r>(-0.6)]<-NA
#   mat_cor_r[mat_cor_r>(-0.6)& mat_cor_r!=1]<-0
#   row.names(mat_cor_r)<-rep(" ",1581)
#   row.names(mat_cor_r)[c(49,283,530,792,1071,1367)]<-c(500,1001,1500,2000,2501,3000)
#   colnames(mat_cor_r)<-rep(" ",1581)
#   colnames(mat_cor_r)[c(49,283,530,792,1071,1367)]<-c(500,1001,1500,2000,2501,3000)
#   #pdf(paste(outpath_heatmap_IRCN_neg0pos0_anno,name_mat,"_IRCN_pheatmap.pdf",sep=""),
#   #    height =6, width = 10)
#   #CairoPNG(file = paste(outpath_heatmap_IRCN_neg6_anno,name_mat,"_neg6_IRCN_pheatmap.png",sep=""), width=8.5, height=6, bg="white")
#   device="png"
#   Cairo(file = paste(outpath_heatmap_IRCN_neg6_anno,name_dmat,"_neg6_IRCN_pheatmap", ".",device , sep=""), 
#         unit="in", dpi=300, width=8.5, height=6, type=device, bg="white")
#   #mat_cor_r = mat_cor_r[rowSums(!is.na(mat_cor_r))!=0, colSums(!is.na(mat_cor_r))!=0]
#   p<-Heatmap(mat_cor_r,
#              #mat_cor_r[1:50,1:50],
#              na_col = "#FFFFFFFF",
#              name = "PCC",
#              col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
#              cluster_columns=F,
#              cluster_rows = F,
#              #column_title = gsub("CC124__","",name_mat), 
#              #column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#              #column_title_side = "top",
#              show_row_names = T,
#              show_column_names = T,
#              #border = T,
#              row_names_side = "left",
#              column_names_side = "bottom",
#              #left_annotation = ha,
#              #row_names_gp = gpar(fontsize = 2),
#              #bottom_annotation = ha,
#              #row_names_gp = gpar(fontsize = 5),
#              #column_names_gp = gpar(fontsize = 5),
#              column_names_rot = 0,
#              gap = unit(0.05, "cm"),
#              #border_gp=gpar(lwd=0.1),
#              #left_annotation = T,
#              left_annotation = rowAnnotation(Mean = anno_lines(as.numeric(Mean_SCRS),
#                                                                
#                                                                axis_param=list(side="top",labels_rot=0),
#                                                                width = unit(1, "cm")),
#                                              #HI = anno_lines(as.numeric(HI_SCRS)),
#                                              #Degree = anno_lines(as.numeric(Degree_SCRS[name_mat,])),
#                                              HI = anno_barplot(
#                                                as.numeric(HI_SCRS), 
#                                                bar_width = 0.0001, 
#                                                gp = gpar(col = "Grey20", fill = "Grey20"), 
#                                                border = T,
#                                                border_gp=gpar(lwd=0.0001,col="Grey20"),
#                                                axis_param=list(side="top",labels_rot=0),
#                                                width = unit(1, "cm")),
#                                              Degree = anno_barplot(
#                                                as.numeric(Degree_SCRS[name_mat,]), 
#                                                bar_width = 0.0001, 
#                                                gp = gpar(col = "#66CCFF", fill = "#66CCFF"), 
#                                                
#                                                axis_param=list(side="top",labels_rot=0),
#                                                border_gp=gpar(lwd=0.0001,col="Grey20"),
#                                                width = unit(2, "cm")),
#                                              annotation_name_side="top",
#                                              #annotation_name_rot=0,
#                                              
#                                              gap = unit(0.2, "cm")),
#              show_heatmap_legend = F
#   )
#   draw(p)
#   decorate_heatmap_body("PCC", {
#     i = which(colnames(mat_cor_p) == "B394")
#     j = which(colnames(mat_cor_p) == "B3341")
#     x = i/ncol(mat_cor_p)
#     y = j/ncol(mat_cor_p)
#     #grid.lines(c(x, x), c(x, y), gp = gpar(lwd = 1, lty = 1,col="black"))
#     grid.polyline(#x=c(y, x), y=c(x, y),  
#       x=c(x,x,y,x),y=c(x,y,x,x),
#       gp=gpar(col="black", lwd=1))
#     #line(c(y, 0), c(0, y), gp = gpar(lwd = 1, lty = 1,col="Red"))
#     #grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
#   })
#   
#   #print(p)
#   dev.off()
#   
# }


#-------------------------------
# Graph cluster based on connectedness (correlation and p-value)
#-------------------------------
# 1. cor_mats_rpvalue format change
v_cor_mats_r_df<-v_cor_mats_rpvalue[[1]][,1:2]
v_cor_mats_pvalue_df<-v_cor_mats_rpvalue[[1]][,1:2]
for(name in names(v_cor_mats_rpvalue)){
  #name<-"000h"
  v_cor_mats_r_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"value"])
  v_cor_mats_pvalue_df_tem<-data.frame(v_cor_mats_rpvalue[[name]][,"p_value"])
  colnames(v_cor_mats_r_df_tem)<-colnames(v_cor_mats_pvalue_df_tem)<-name
  v_cor_mats_r_df<-cbind(v_cor_mats_r_df,v_cor_mats_r_df_tem)
  v_cor_mats_pvalue_df<-cbind(v_cor_mats_pvalue_df,v_cor_mats_pvalue_df_tem)
}
#v_cor_mats_pvalue_df[1:5,15:16]
write.csv(file=paste(outpath,'global_connectness_v_cor_mats_r_df.csv',sep=""),as.data.frame(v_cor_mats_r_df), 
          quote=F, row.names=T)
write.csv(file=paste(outpath,'global_connectness_v_cor_mats_pvalue_df.csv',sep=""),as.data.frame(v_cor_mats_pvalue_df), 
          quote=F, row.names=T)
# 2. v_cor_mats_df ramanomes distance (connectedness all & p.value no cutoff)
require("proxy")
dist_cor_mats_r<-as.matrix(dist(t(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)]), method="Euclidean")) #Euclidean dist
#dist_cor_mats<-as.matrix(dist(t(v_cor_mats_df[1:10,7:22]), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'global_connectness_euclidean_dist_cor_mats_r_hclust.csv',sep=""),as.data.frame(dist_cor_mats_r), 
          quote=F, row.names=T)
#https://www.jianshu.com/p/d6b361c479d1
mtrx2cols_1col = function(m1,val1){
  lt = lower.tri(m1)  
  res = data.frame(row = row(m1,as.factor = T)[lt],  
                   col = col(m1,as.factor = T)[lt],  
                   val1 = m1[lt]) 
  names(res)[3] = c(val1) 
  return(res)
}
dist_cor_mats_r_df<-mtrx2cols_1col(dist_cor_mats_r,"Euclidean_dist")
write.csv(file=paste(outpath,'global_connectness_euclidean_dist_cor_r_df_hclust.csv',sep=""),as.data.frame(dist_cor_mats_r_df), 
          quote=F, row.names=T)
# Graph cluster based on dist_cor_mats
CairoPDF(file = paste(outpath, "global_connectness_euclidean_dist_cor_r_df_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
#pdf(file = paste(outpath,'global_connectness_euclidean_dist_cor_df_hclust.pdf',sep=""), width=3+nlevels(Group)*1.5, height=70);
plot(hclust(as.dist(dist_cor_mats_r), "ward.D"))
dev.off()

# 3.v_cor_mats_rpvalue_df ramanomes distance(connectedness all & p.value<=0.05)
v_cor_mats_rpvalue_df<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
#v_cor_mats_pvalue_df_new<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
#v_cor_mats_rpvalue_df[which(v_cor_mats_pvalue_df_new>0.05)]<-0 # p.value>0.05
#v_cor_mats_pvalue_df_new_tem<-v_cor_mats_pvalue_df_new[,15:16]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_pvalue_df_partial.csv',sep=""),as.data.frame(v_cor_mats_pvalue_df_new_tem), 
#          quote=F, row.names=T)
#v_cor_mats_rpvalue_df_tem<-v_cor_mats_rpvalue_df[,15:16]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_rpvalue_df_partial.csv',sep=""),as.data.frame(v_cor_mats_rpvalue_df_tem), 
#          quote=F, row.names=T)
dist_cor_mats_rpvalue<-as.matrix(dist(t(v_cor_mats_rpvalue_df), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'global_connectness_euclidean_dist_cor_mats_rpvalue_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "global_connectness_euclidean_dist_cor_rpvalue_df_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
#pdf(file = paste(outpath,'global_connectness_euclidean_dist_cor_df_neg0.6_hclust.pdf',sep=""), width=3+nlevels(Group)*1.5, height=70);
plot(hclust(as.dist(dist_cor_mats_rpvalue), "ward.D"))
dev.off()


# 4.v_cor_mats_rpvalue_df ramanomes distance(connectedness<-0.6 & p.value<=0.05)
v_cor_mats_rpvalue_df_neg6<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_pvalue_df_neg6<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
#v_cor_mats_rpvalue_df_neg6[which(v_cor_mats_rpvalue_df_neg6>(-0.6) | v_cor_mats_pvalue_df_neg6>(0.05))]<-0 #rho>-0.6 or p.value>0.05
v_cor_mats_rpvalue_df_neg6[which(v_cor_mats_rpvalue_df_neg6>(-0.6))]<-0 #rho>-0.6 
v_cor_mats_rpvalue_df_neg6[which(v_cor_mats_pvalue_df_neg6>(0.05))]<-0 #p.value>0.05
#v_cor_mats_pvalue_df_neg6_tem<-v_cor_mats_pvalue_df_neg6[,15:16]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_pvalue_df_neg6_partial.csv',sep=""),as.data.frame(v_cor_mats_pvalue_df_neg6_tem), 
#          quote=F, row.names=T)
#v_cor_mats_rpvalue_df_neg6_tem<-v_cor_mats_rpvalue_df_neg6[,15:16]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_rpvalue_df_neg6_partial.csv',sep=""),as.data.frame(v_cor_mats_rpvalue_df_neg6_tem), 
#          quote=F, row.names=T)
dist_cor_mats_rpvalue_neg6<-as.matrix(dist(t(v_cor_mats_rpvalue_df_neg6), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'global_connectness_euclidean_dist_cor_mats_rpvalue_neg6_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue_neg6), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "global_connectness_euclidean_dist_cor_rpvalue_df_neg0.6_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
#pdf(file = paste(outpath,'global_connectness_euclidean_dist_cor_df_neg0.6_hclust.pdf',sep=""), width=3+nlevels(Group)*1.5, height=70);
plot(hclust(as.dist(dist_cor_mats_rpvalue_neg6), "ward.D"))
dev.off()

# 5.v_cor_mats_rpvalue_df ramanomes distance(connectedness<-0.8 & p.value<=0.05)
v_cor_mats_rpvalue_df_neg8<-as.matrix(v_cor_mats_r_df[,3:ncol(v_cor_mats_r_df)])
v_cor_mats_pvalue_df_neg8<-as.matrix(v_cor_mats_pvalue_df[,3:ncol(v_cor_mats_pvalue_df)])
#v_cor_mats_rpvalue_df_neg8[which(v_cor_mats_rpvalue_df_neg8>(-0.8) | v_cor_mats_pvalue_df_neg8>0.05)]<-0 #rho>-0.8 or p.value>0.05
v_cor_mats_rpvalue_df_neg8[which(v_cor_mats_rpvalue_df_neg8>(-0.8))]<-0 #rho>-0.6 
v_cor_mats_rpvalue_df_neg8[which(v_cor_mats_pvalue_df_neg8>(0.05))]<-0 #p.value>0.05
#v_cor_mats_pvalue_df_neg8_tem<-v_cor_mats_pvalue_df_neg8[,15:16]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_pvalue_df_neg8_partial.csv',sep=""),as.data.frame(v_cor_mats_pvalue_df_neg8_tem), 
#          quote=F, row.names=T)
#v_cor_mats_rpvalue_df_neg8_tem<-v_cor_mats_rpvalue_df_neg8[,1:20]
#write.csv(file=paste(outpath,'global_connectness_v_cor_mats_rpvalue_df_neg8_partial.csv',sep=""),as.data.frame(v_cor_mats_rpvalue_df_neg8_tem), 
#          quote=F, row.names=T)
dist_cor_mats_rpvalue_neg8<-as.matrix(dist(t(v_cor_mats_rpvalue_df_neg8), method="Euclidean")) #Euclidean dist
write.csv(file=paste(outpath,'global_connectness_euclidean_dist_cor_mats_rpvalue_neg8_hclust.csv',sep=""),as.data.frame(dist_cor_mats_rpvalue_neg8), 
          quote=F, row.names=T)
device="pdf"
CairoPDF(file = paste(outpath, "global_connectness_euclidean_dist_cor_rpvalue_df_neg0.8_hclust.pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
#pdf(file = paste(outpath,'global_connectness_euclidean_dist_cor_df_neg0.6_hclust.pdf',sep=""), width=3+nlevels(Group)*1.5, height=70);
plot(hclust(as.dist(dist_cor_mats_rpvalue_neg8), "ward.D"))
dev.off()

#-------------------------------
# Graph cluster based on mean spectra
#-------------------------------
mean_spectra<-aggregate(mat,list(Group), mean)
rownames(mean_spectra)<-mean_spectra[,1]
mean_spectra<-mean_spectra[,-1]
write.csv(data.frame(mean_spectra),paste(outpath,"meanspectra_mats_",matrixfile,".csv",sep=""),row.names=T,quote=F)
#write.table(mean_spectra,paste(outpath,"meanspectra_mats_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
mean_spectra_mat<-as.matrix(mean_spectra)
# 1. jsd distance
jsd_meanspectra<-JSD(mean_spectra_mat)
jsd_meanspectra[is.na(jsd_meanspectra)]<-0
#write.table(jsd_meanspectra,paste(outpath,"jsd_meanspectra_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
write.csv(data.frame(jsd_meanspectra),paste(outpath,"jsd_meanspectra_",matrixfile,".csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "global_jsd_meanspectra_hclust.", category, ".pdf" , sep=""),width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(jsd_meanspectra), "ward.D"))
dev.off()
# 2. cosine distance
require("proxy")
cosine_meanspectra<-as.matrix(proxy::dist(mean_spectra_mat, method="cosine")) #cosine dist
cosine_meanspectra[is.na(cosine_meanspectra)]<-0
#write.table(data.frame(cosine_meanspectra),paste(outpath,"cosine_meanspectra_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
write.csv(data.frame(cosine_meanspectra),paste(outpath,"cosine_meanspectra_",matrixfile,".csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "global_cosine_meanspectra_hclust.", category, ".pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(cosine_meanspectra), "ward.D"))
dev.off()
# 3. euclidean distance
require("proxy")
euclidean_meanspectra<-as.matrix(proxy::dist(mean_spectra_mat, method="Euclidean")) # euclidean dist
euclidean_meanspectra[is.na(euclidean_meanspectra)]<-0
#write.table(data.frame(euclidean_meanspectra),paste(outpath,"euclidean_meanspectra_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
write.csv(data.frame(euclidean_meanspectra),paste(outpath,"euclidean_meanspectra_",matrixfile,".csv",sep=""),row.names=T,quote=F)
CairoPDF(file = paste(outpath, "global_euclidean_meanspectra_hclust.", category, ".pdf", sep=""), width=3+nlevels(Group)*1.5, height=70, bg="white");
plot(hclust(as.dist(euclidean_meanspectra), "ward.D"))
dev.off()
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------
#Bands marker selection by degree
#-------------------------------
#E:\RWAS\20201116_NetworkAnalysis_CC124\NetworkAnalysis_Bands_Degree_ggbarplot_20201117.r
#bands_Cluster<-read.csv(paste(input,
#                              "global_network_Cluster_membership_by_Celltype__Timepoint.xls",sep = ""),
#                      sep = "\t" )
#bands_Degree<-read.csv(paste(input,
#                              "global_network_Degree_by_Celltype__Timepoint.xls",sep = ""),
#                        sep = "\t" )
#bands_Cluster<-Cluster_walktrap_membership_df
#bands_Degree<-Degree_df
# ggbarplot
#output_ggdotchart<-paste("./","ggdotchart/",sep = "")
#dir.create(output_ggdotchart)
#for (name_col in colnames(bands_Degree)[5:20]) {
  #dfm<-data.frame(Wave_num=bands_Degree$Wave_num,Degree=bands_Degree$CC124__012h,Cluster=bands_Cluster$CC124__012h)
  #name_col<-"CC124__000h"
 # dfm<-data.frame(Wave_num=bands_Degree$Wave_num,
  #                Degree=bands_Degree[,name_col],
  #                Cluster=bands_Cluster[,name_col])
  #order_degree<-order(dfm$Degree,decreasing = T)
  #dfm<-dfm[order_degree,]
  
 # Freq_cluster<-data.frame(table(dfm$Cluster))
 # order_cluster<-order(Freq_cluster$Freq,decreasing = T)
 # Freq_cluster<-Freq_cluster[order_cluster,]
 # dfm$Cluster_size<-ifelse(dfm$Cluster==Freq_cluster$Var1[1]& dfm$Degree!=0,"Cluster_1st",
  #                         ifelse(dfm$Cluster==Freq_cluster$Var1[2]& dfm$Degree!=0,"Cluster_2nd",
  #                                ifelse(dfm$Cluster==Freq_cluster$Var1[3]& dfm$Degree!=0,"Cluster_3rd","NA")))
  
 # dfm_sub_1<-dfm[which(dfm$Cluster_size=="Cluster_1st"),]
 # dfm_sub_1<-if(nrow(dfm_sub_1)<10){
 #   dfm_sub_1
 # } else{
 #   dfm_sub_1[1:10,]
 # }
  
 # dfm_sub_2<-dfm[which(dfm$Cluster_size=="Cluster_2nd"),]
 # dfm_sub_2<-if(nrow(dfm_sub_2)<10){
 #   dfm_sub_2
 # } else{
 #   dfm_sub_2[1:10,]
 # }
  
#  dfm_sub_3<-dfm[which(dfm$Cluster_size=="Cluster_3rd"),]
#  dfm_sub_3<-if(nrow(dfm_sub_3)<10){
#    dfm_sub_3
#  } else{
#    dfm_sub_3[1:10,]
#  }
#  dfm_sub<-rbind(dfm_sub_1,dfm_sub_2,dfm_sub_3)
  
  #dfm_sub<-dfm_sub[order(dfm_sub$Cluster_size,decreasing = F),]
#  dfm_sub$Cluster_size<-factor(dfm_sub$Cluster_size,
 #                              levels = c("Cluster_3rd","Cluster_2nd","Cluster_1st"))
  #cluster_cols <- c("#00AFBB", "#E7B800", "#FC4E07")
  
#  #dfm_sub1<-dfm_sub[nrow(dfm_sub):1,]
#  require(ggpubr)
#  library(ggpubr)
#  dfm_sub$Wave_num<-as.character(dfm_sub$Wave_num)
#  p<- ggdotchart(dfm_sub, x = "Wave_num", y = "Degree",
#               color = "Cluster_size",                                # Color by groups
#               palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
#               sorting = "descending",                       # Sort value in descending order
#               add = "segments",                             # Add segments from y = 0 to dots
#               rotate = TRUE,                                # Rotate vertically
#               group = "Cluster_size",                                # Order by groups
#               dot.size = 4,                                 # Large dot size
#               label = round(dfm_sub$Degree),                        # Add mpg values as dot labels
#               font.label = list(color = "white", size =4, 
#                                 vjust = 0.5),               # Adjust label parameters
#               ggtheme = theme_pubr()                      # ggplot2 theme#    )
#    )
  #ggsave(paste(output,"CC124__018h","_Bands_Degree_ggbarplot.pdf",sep = ""),
  #       p,width = 4,height = 0.2*nrow(dfm_sub)+0.1)
#  ggsave(paste(output_ggdotchart,name_col,"_Bands_Degree_ggbarplot.pdf",sep = ""),
#         p,width = 4,height = 0.22*nrow(dfm_sub)+0.1)
#}
#----------------------------------------------------------------------------------------------------

#-------------------------------
#PCA based on MP/MI/MC
#-------------------------------
#E:\RWAS\20201116_NetworkAnalysis_CC124\NetworkAnalysis_PCA_MPMIMC_20201117.r
#--------------
# 1. MP PCA mean_spectra
#--------------
mat<-mean_spectra
tem<-mat[1:2,1:5]
#mat<-mat[which(row.names(mat)%in%sub_title),]
MP.pca <- prcomp(mat, scale = F)
#--------------
# 2. MC PCA (connectedness<-0.6 & p.value<=0.05)
#--------------
mat<-t(v_cor_mats_rpvalue_df_neg6)
tem<-mat[1:2,1:5]
#mat<-mat[which(row.names(mat)%in%sub_title),]
MC.pca <- prcomp(mat, scale = F)
#-------------
# 3. MI PCA (connectedness all & p.value<=0.05)
#--------------
mat<-t(v_cor_mats_rpvalue_df)
tem<-mat[1:2,1:5]
tem<-mat[1:2,ncol(mat)]
tem<-mat[1:2,(ncol(mat)-2):ncol(mat)]
#mat<-mat[which(row.names(mat)%in%sub_title),]
MI.pca <- prcomp(mat, scale = F)
#--------------
# 4. PCA.RData
#--------------
save(MP.pca, MC.pca, MI.pca, file = paste(outpath,"PCA_MP-MC-MI.RData",sep=""))
#save(MP.pca,file = paste(outpath,"PCA_MP.RData",sep=""))
#load(paste(input,"PCA_MP-MC-MI.RData.RData",sep = ""))

# Desktopï¼PCA
##########################
#file.path<-"E:/RWAS/20201116_NetworkAnalysis_CC124/"
#setwd(file.path)
#input <-paste("E:/RWAS/Network_20200514/Timepoint_neg0.6/",sep="")
output_PCA <- paste(outpath,"./Results_PCA_MPMIMC_20201117/",sep="")
dir.create(output_PCA)

#load(paste(input,"NoPos-Neg_0.6_PCA_MP-MC-MI.RData",sep = ""))
#Group_abbr<-read.table(paste("Name_abbr.txt",sep = ""),header = T, sep = "\t")
library(factoextra)
#--------------
# 1. MP PCA 
#load(paste(input,"PCA_MP.RData",sep = ""))
#Group_cluster<-read.csv(paste("E:/RWAS/Network_20200525/Results_20201117/",
#                              "Cluster_6_meanspectra.csv",sep = ""))
A<-fviz_eig(MP.pca, addlabels = TRUE, ylim = c(0, 50))
#A$eig[1:10]
#b<-MC.pca$x
#Group_abbr$pca_x_rowname<-row.names(MP.pca$x)
#row.names(MP.pca$x)<-Group_abbr$Group_abbr
#groups<-as.factor(Group_cluster$Cluster)
pdf(paste(output_PCA,"PCA_MP.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MP.pca,
             xlab="PC1 (41.84%)",
             ylab="PC2 (28.55%)",
#             col.ind = groups,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
#--------------
# 2. MC PCA 
#Group_cluster<-read.csv(paste("E:/RWAS/Network_20200525/Results_20201117/",
#                              "Cluster_6_dist_cor_mats_rpvalue_neg6.csv",sep = ""))
A<-fviz_eig(MC.pca, addlabels = TRUE, ylim = c(0, 50))
#A$eig[1:10]
#b<-MC.pca$x
#Group_abbr$pca_x_rowname<-row.names(MP.pca$x)
#row.names(MC.pca$x)<-Group_abbr$Group_abbr
#groups<-as.factor(Group_cluster$Cluster)
pdf(paste(output_PCA,"PCA_MC.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MC.pca,
             xlab="PC1 (22.07%)",
             ylab="PC2 (12.26%)",
#             col.ind = groups,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
#--------------
# 3. MI PCA 
#Group_cluster<-read.csv(paste("E:/RWAS/Network_20200525/Results_20201117/",
#                              "Cluster_6_dist_cor_mats_rpvalue0.csv",sep = ""))
A<-fviz_eig(MI.pca, addlabels = TRUE, ylim = c(0, 50))
#A$eig[1:10]
#b<-MC.pca$x
#Group_abbr$pca_x_rowname<-row.names(MP.pca$x)
#row.names(MI.pca$x)<-Group_abbr$Group_abbr
#groups<-as.factor(Group_cluster$Cluster)
pdf(paste(output_PCA,"PCA_MI.pdf",sep=""),height = 8,width = 8)
fviz_pca_ind(MI.pca,
             xlab="PC1 (25.38%)",
             ylab="PC2 (15.57%)",
#             col.ind = groups,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
#------------------------------------------------------


#-------------------------------
# Local IRCNs by chordDiagram
# Generate local networks by selected bands
#-------------------------------
#Peak_sub<-c("B907","B959","B1005","B1129",
#         "B1156",  "B1186", "B1260","B1489")
#Peak_sub<-paste("B",local_bands_ann$Wave_num,sep = "")
Peak_sub<-as.character(local_bands_ann$Wave_num)
bands<-Peak_sub
#aaa<-cor_mats_rpvalue[["Yeast__day 0"]]$r
#which(bands%in%colnames(aaa))
Peaks_fixed_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, 
                                   function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
names(Peaks_fixed_cor_mats_rpvalue)<-names(cor_mats_rpvalue)
save(Peaks_fixed_cor_mats_rpvalue,
     file=paste(outpath,
     #file=paste("/mnt/data6/heyh/20180831_Analysis_RWAS/Group_neg0.6/",
                "Peaks_fixed_mannulselected_cor_mats_rpvalue.RData",sep=""))

#------------------------------------------------------------------

#------------------------------------------------------------------
if(!is.null(LocalBandsAnnfile) & plot_local_IRCN){ 
  
  #Peak_sub<-paste("B",local_bands_ann$Wave_num,sep = "")
  Peak_sub<-local_bands_ann$Wave_num
  #-------------------------------
  # Generate local networks by selected bands
  #-------------------------------
  bands<-as.character(Peak_sub)
  bands_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, 
                                   function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
  require(Cairo)
  local_bands_ann_new<-local_bands_ann
local_bands_ann_new$Wave_num<-as.character(local_bands_ann_new$Wave_num)
local_bands_ann_new$Wave_num<-as.character(round(as.numeric(gsub("B","",as.character(local_bands_ann_new$Wave_num))),0))

  #---------------------------------------------------
  # png plot_local_chorddiagram
  #---------------------------------------------------
  # 1.# Pos_Edge=F, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_neg6_png_rpvalue <-paste(outpath_png,"local_neg6_png_rpvalue/",sep = "")
  dir.create(outpath_local_neg6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
	#group_name<-"CC124__002h"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    #CairoPNG(file = paste(outpath_local_neg6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
    Cairo(file = paste(outpath_local_neg6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0.6, local_bands_ann_new)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 2.# Pos_Edge=T, Neg_Edge=F, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_pos6_png_rpvalue <- paste(outpath_png,"local_pos6_png_rpvalue/",sep = "")
  dir.create(outpath_local_pos6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #CairoPNG(file = paste(outpath_local_pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""),width=3, height=3, bg="white")
	Cairo(file = paste(outpath_local_pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0.6, local_bands_ann_new)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 3.# Pos_Edge=T, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_neg6pos6_png_rpvalue <- paste(outpath_png,"local_neg6pos6_png_rpvalue/",sep = "")
  dir.create(outpath_local_neg6pos6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #CairoPNG(file = paste(outpath_local_neg6pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
	Cairo(file = paste(outpath_local_neg6pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0.6, local_bands_ann_new)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 4.# Pos_Edge=F, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_neg0_png_rpvalue <-  paste(outpath_png,"local_neg0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0_png_rpvalue/"
  dir.create(outpath_local_neg0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #CairoPNG(file = paste(outpath_local_neg0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png" , sep=""), width=3, height=3, bg="white")
	Cairo(file = paste(outpath_local_neg0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0, local_bands_ann_new)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 5.# Pos_Edge=T, Neg_Edge=F, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_pos0_png_rpvalue <- paste(outpath_png,"local_pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_pos0_png_rpvalue/"
  dir.create(outpath_local_pos0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #CairoPNG(file = paste(outpath_local_pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, sep=""), width=3, height=3, bg="white")
	Cairo(file = paste(outpath_local_pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0, local_bands_ann_new)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 6.# Pos_Edge=T, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_neg0pos0_png_rpvalue <- paste(outpath_png,"local_neg0pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0pos0_png_rpvalue/"
  dir.create(outpath_local_neg0pos0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #CairoPNG(file = paste(outpath_local_neg0pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
	Cairo(file = paste(outpath_local_neg0pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".",device , sep=""), 
        unit="in", dpi=300, width=3, height=3, type=device, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0, local_bands_ann_new)
    dev.off()  
  }
  #----------------------------------------------------------------------------------------------
}
#---------------------------------------------------------------------------------------------------------------------


#--------------------------------
# save data
#--------------------------------
#save(cor_mats,v_cor_mats,v_cor_mats_df,v_cor_mats_trimmed,g_array,g,Degree,mean_spectra,dist_cor_mats,
#     file = paste(outpath,"mydata_maindataset_by_", category,".RData",sep=""))
save.image(file = paste(outpath,"mydata_by_", category,".RData",sep="")) #save all data into one .RData file


#-------------------------------
# Local IRCNs by chordDiagram
# Desk top
#-------------------------------
#---------------------
# loading the dataset
#---------------------
#load(paste("/mnt/data6/heyh/IRCN_CC124_20200624/Timepoint_neg0.6/",
#           "mydata_by_Celltype__Timepoint.RData",sep=""))
# load(paste("E:/RWAS/IRCN_CC124_DU_20201022/Timepoint_neg0.6/",
#            "Peaks_fixed_mannulselected_cor_mats_rpvalue.RData",sep=""))

bands_cor_mats_rpvalue<-Peaks_fixed_cor_mats_rpvalue

# # global_bands_annotation é·æ³ç?(+-3)
# Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
#                                           "Global_bands_annotation_new.txt",sep = ""),
#                                     header = T,sep="\t")
# i<-4
# #i<-79
# while (i<=(nrow(Global_bands_annotation)-3)){
#   if(Global_bands_annotation$Group_new[i]!="Unknown"){
#     Global_bands_annotation$Group_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Group_new[i]),6)
#     Global_bands_annotation$Assignment_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Assignment_new[i]),6)
#     i<-i+4
#   }else{
#       i<-i+1
#     }
# }
# write.csv(data.frame(Global_bands_annotation),
#           paste("E:/RWAS/IRCN_CC124_DU_20201022/","Global_bands_annotation_expanded.csv",sep = ""),
#           quote = F,row.names = F)
# 
# 
# 
# Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
#                                           "Global_bands_annotation.txt",sep = ""),
#                                     header = T,sep="\t")
# i<-4
# #i<-79
# Global_bands_annotation$Group_new<-Global_bands_annotation$Group
# Global_bands_annotation$Assignment<-Global_bands_annotation$Assignment_new
# while (i<=(nrow(Global_bands_annotation)-3)){
#   if(Global_bands_annotation$Group[i]!="Unknown"){
#     Global_bands_annotation$Group_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Group_new[i]),6)
#     Global_bands_annotation$Assignment_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Assignment_new[i]),6)
#     i<-i+4
#   }else{
#     i<-i+1
#   }
# }

# Libraries
library(dplyr)
library(tidyverse)
library(circlize)
options(knitr.table.format = "html")
library(viridis)
library(igraph)

library(colormap)
library(hrbrthemes)
library(kableExtra)
library(ggraph)#graphlayouts #tidygraph


################
#1.Neg 0 & p.value<0.05 all peaks
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath_Arc <- paste(outpath,"1.Arcdiagram_Peaks_Neg0/",sep="")
#outpath <- "E:/RWAS/IRCN_CC124_DU_20201022/Timepoint_neg0.6/1.Arcdiagram_Peaks_Neg0+/"
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  #Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(0)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  #----------------------
  #add group
  # local_bands_ann<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
  #                                   "Local_bands_annotation_DU.txt",sep = ""),
  #                             header = T,sep="\t")
  # local_bands_ann$Wave_num<-paste("B",local_bands_ann$Wave_num,sep = "")
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
  # #--------------------
  # # 2. Make the graph
  # # mannul grouped
  # #--------------------
  # #order_Peaks<-c("B957","B966","B999","B1007",
  # #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  # #order_Peaks<-c("B1526","B1158","B1007",
  # #               "B1504","B1152","B999",
  # #               "B1193","B1180","B957","B966")
  # order_Peaks<-c("B478","B865","B940","B1049","B1049","B1049","B1049",
  #                "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
  #                "B526",  "B577", "B612","B718")
  # order_Peaks<-c(#"B16581441",
  #                "B1007","B1602",
  #                "B478","B865","B938","B1125",
  #                "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
  #                "B526",  "B577", 
  #                "B612","B718",
  #                "B2911")
  # # Transform the adjacency matrix in a long format
  # connect <- dataUU %>% 
  #   gather(key="to", value="value", -1) %>%
  #   mutate(to = gsub("\\.", " ",to)) %>%
  #   na.omit() 
  # 
  # # Number of connection per person
  # c(as.character(connect$from), as.character(connect$to)) %>%
  #   as.tibble() %>%
  #   group_by(value) %>%
  #   summarize(n=n()) -> coauth
  # 
  # colnames(coauth) <- c("name", "n")
  # bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  # bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  # coauth<-rbind(coauth,bands_more)
  # 
  # #dim(coauth)
  # new_order<-match(order_Peaks,coauth$name)
  # coauth<-coauth[new_order,]
  # 
  # # grouped
  # #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  # #                                    header = T,sep="\t")
  # Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  # new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  # Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  # 
  # coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  # coauth$grp_mannul<-c(#"DU",
  #                      rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
  #                      rep("Lipids",2),"Lipids; Carbohydrates")
  # coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(##"DU",
  #                                                      "Proteins","Starch","TAG", "Memanbrane Lipids",
  #                                                      "Lipids","Lipids; Carbohydrates"))
  # #coauth$grp_mannul<-c(#"DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU","Starch","TAG", "Memanbrane Lipids"))
  # #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,
  # #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  # 
  # #coauth$grp_mannul<-coauth$Group_new
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  # #coauth<-coauth[order(coauth$grp_mannul),]
  # 
  # # keep only this people in edges
  # connect <- connect %>%
  #   filter(from %in% coauth$name) %>%
  #   filter(to %in% coauth$name)
  # 
  # # Create a graph object with igraph
  # mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  # 
  # # prepare a vector of n color in the viridis scale
  # mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  # #mycolor <- sample(mycolor, length(mycolor))
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  # #mycolor<-c("#440154ff", "#5cc863ff")
  # #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  # 
  # #----------------------
  # #no group and assignment
  # p21<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0.4,0), "null"),
  #     panel.spacing=unit(c(0,0,3.4,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  # ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  # 
  # #----------------------
  # #add assignment
  # p22<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
  #        p22,height = 3,width = 3)
  # 
  # #----------------------
  # #add group
  # p23<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
  #        p23,height = 4,width = 2)
  # 
  # 
  }


################
#2.Pos 0 & p.value<0.05 
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath_Arc <- paste(outpath,"2.Arcdiagram_Peaks_Pos0/",sep="")
#outpath <- "E:/RWAS/IRCN_CC124_DU_20201022/Timepoint_neg0.6/1.Arcdiagram_Peaks_Neg0+/"
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat<(0)]<-NA # r<0
  #Corr_mat[which(colnames(Corr_mat)!="B16581441"),which(colnames(Corr_mat)!="B16581441")]<-NA # DU only
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  #----------------------
  #add group
  # local_bands_ann<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
  #                                   "Local_bands_annotation_DU.txt",sep = ""),
  #                             header = T,sep="\t")
  # local_bands_ann$Wave_num<-paste("B",local_bands_ann$Wave_num,sep = "")
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  # 
  # #--------------------
  # # 2. Make the graph
  # # mannul grouped
  # #--------------------
  # #order_Peaks<-c("B957","B966","B999","B1007",
  # #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  # #order_Peaks<-c("B1526","B1158","B1007",
  # #               "B1504","B1152","B999",
  # #               "B1193","B1180","B957","B966")
  # order_Peaks<-c(#"B16581441",
  #                "B1007","B1602",
  #                "B478","B865","B938","B1125",
  #                "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
  #                "B526",  "B577", 
  #                "B612","B718",
  #                "B2911")
  # # Transform the adjacency matrix in a long format
  # connect <- dataUU %>% 
  #   gather(key="to", value="value", -1) %>%
  #   mutate(to = gsub("\\.", " ",to)) %>%
  #   na.omit() 
  # 
  # # Number of connection per person
  # c(as.character(connect$from), as.character(connect$to)) %>%
  #   as.tibble() %>%
  #   group_by(value) %>%
  #   summarize(n=n()) -> coauth
  # 
  # colnames(coauth) <- c("name", "n")
  # bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  # bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  # coauth<-rbind(coauth,bands_more)
  # 
  # #dim(coauth)
  # new_order<-match(order_Peaks,coauth$name)
  # coauth<-coauth[new_order,]
  # 
  # # grouped
  # #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  # #                                    header = T,sep="\t")
  # Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  # new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  # Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  # 
  # coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  # coauth$grp_mannul<-c(#"DU",
  #   rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
  #                      rep("Lipids",2),"Lipids; Carbohydrates")
  # coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
  #                                                      "Proteins","Starch","TAG", "Memanbrane Lipids",
  #                                                      "Lipids","Lipids; Carbohydrates"))
  # #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  # #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,
  # #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  # 
  # #coauth$grp_mannul<-coauth$Group_new
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  # #coauth<-coauth[order(coauth$grp_mannul),]
  # 
  # # keep only this people in edges
  # connect <- connect %>%
  #   filter(from %in% coauth$name) %>%
  #   filter(to %in% coauth$name)
  # 
  # # Create a graph object with igraph
  # mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  # 
  # # prepare a vector of n color in the viridis scale
  # mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  # #mycolor <- sample(mycolor, length(mycolor))
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  # #mycolor<-c("#440154ff", "#5cc863ff")
  # #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  # 
  # #----------------------
  # #no group and assignment
  # p21<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0.4,0), "null"),
  #     panel.spacing=unit(c(0,0,3.4,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  # ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  # 
  # #----------------------
  # #add assignment
  # p22<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
  #        p22,height = 3,width = 3)
  # 
  # #----------------------
  # #add group
  # p23<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
  #        p23,height = 4,width = 2)
  # 
  
}

################
#3.Neg 6 & p.value<0.05 all peaks
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"

outpath_Arc <- paste(outpath,"3.Arcdiagram_Peaks_Neg6/",sep="")
#outpath <- "E:/RWAS/IRCN_CC124_DU_20201022/Timepoint_neg0.6/1.Arcdiagram_Peaks_Neg0+/"
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  #Corr_mat[Corr_mat>(0)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  #----------------------
  #add group
  # local_bands_ann<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
  #                                   "Local_bands_annotation_DU.txt",sep = ""),
  #                             header = T,sep="\t")
  # local_bands_ann$Wave_num<-paste("B",local_bands_ann$Wave_num,sep = "")
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
  # #--------------------
  # # 2. Make the graph
  # # mannul grouped
  # #--------------------
  # #order_Peaks<-c("B957","B966","B999","B1007",
  # #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  # #order_Peaks<-c("B1526","B1158","B1007",
  # #               "B1504","B1152","B999",
  # #               "B1193","B1180","B957","B966")
  # order_Peaks<-c(#"B16581441",
  #                "B1007","B1602",
  #                "B478","B865","B938","B1125",
  #                "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
  #                "B526",  "B577", 
  #                "B612","B718",
  #                "B2911")
  # # Transform the adjacency matrix in a long format
  # connect <- dataUU %>% 
  #   gather(key="to", value="value", -1) %>%
  #   mutate(to = gsub("\\.", " ",to)) %>%
  #   na.omit() 
  # 
  # # Number of connection per person
  # c(as.character(connect$from), as.character(connect$to)) %>%
  #   as.tibble() %>%
  #   group_by(value) %>%
  #   summarize(n=n()) -> coauth
  # 
  # colnames(coauth) <- c("name", "n")
  # bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  # bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  # coauth<-rbind(coauth,bands_more)
  # 
  # #dim(coauth)
  # new_order<-match(order_Peaks,coauth$name)
  # coauth<-coauth[new_order,]
  # 
  # # grouped
  # #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  # #                                    header = T,sep="\t")
  # Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  # new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  # Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  # 
  # coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  # coauth$grp_mannul<-c(#"DU",
  #                      rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
  #                      rep("Lipids",2),"Lipids; Carbohydrates")
  # coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
  #                                                      "Proteins","Starch","TAG", "Memanbrane Lipids",
  #                                                      "Lipids","Lipids; Carbohydrates"))
  # #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  # #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,
  # #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  # 
  # #coauth$grp_mannul<-coauth$Group_new
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  # #coauth<-coauth[order(coauth$grp_mannul),]
  # 
  # # keep only this people in edges
  # connect <- connect %>%
  #   filter(from %in% coauth$name) %>%
  #   filter(to %in% coauth$name)
  # 
  # # Create a graph object with igraph
  # mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  # 
  # # prepare a vector of n color in the viridis scale
  # mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  # #mycolor <- sample(mycolor, length(mycolor))
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  # #mycolor<-c("#440154ff", "#5cc863ff")
  # #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  # 
  # #----------------------
  # #no group and assignment
  # p21<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0.4,0), "null"),
  #     panel.spacing=unit(c(0,0,3.4,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  # ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  # 
  # #----------------------
  # #add assignment
  # p22<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
  #        p22,height = 3,width = 3)
  # 
  # #----------------------
  # #add group
  # p23<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
  #        p23,height = 4,width = 2)
  # 
  # 
}



################
#4.Neg 0.6 & p.value<0.05 DU only
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath_Arc <- paste(outpath,"4.Arcdiagram_Peaks_Pos6/",sep="")
#outpath <- "E:/RWAS/IRCN_CC124_DU_20201022/Timepoint_neg0.6/1.Arcdiagram_Peaks_Neg0+/"
dir.create(outpath_Arc)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # r<0
  #Corr_mat[which(colnames(Corr_mat)!="B16581441"),which(colnames(Corr_mat)!="B16581441")]<-NA # DU only
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath_Arc,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  #----------------------
  #add group
  # local_bands_ann<-read.table(paste("E:/RWAS/IRCN_CC124_DU_20201022/",
  #                                   "Local_bands_annotation_DU.txt",sep = ""),
  #                             header = T,sep="\t")
  # local_bands_ann$Wave_num<-paste("B",local_bands_ann$Wave_num,sep = "")
  coauth$grp_mannul<-local_bands_ann[which(local_bands_ann$Wave_num==coauth$name),"Group"]
  coauth$grp_mannul<-factor(coauth$grp_mannul)
  coauth<-coauth[order(coauth$grp_mannul),]
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#FF6666", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], 
                    fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath_Arc,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
  # #--------------------
  # # 2. Make the graph
  # # mannul grouped
  # #--------------------
  # #order_Peaks<-c("B957","B966","B999","B1007",
  # #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  # #order_Peaks<-c("B1526","B1158","B1007",
  # #               "B1504","B1152","B999",
  # #               "B1193","B1180","B957","B966")
  # order_Peaks<-c(#"B16581441",
  #                "B1007","B1602",
  #                "B478","B865","B938","B1125",
  #                "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
  #                "B526",  "B577", 
  #                "B612","B718",
  #                "B2911")
  # # Transform the adjacency matrix in a long format
  # connect <- dataUU %>% 
  #   gather(key="to", value="value", -1) %>%
  #   mutate(to = gsub("\\.", " ",to)) %>%
  #   na.omit() 
  # 
  # # Number of connection per person
  # c(as.character(connect$from), as.character(connect$to)) %>%
  #   as.tibble() %>%
  #   group_by(value) %>%
  #   summarize(n=n()) -> coauth
  # 
  # colnames(coauth) <- c("name", "n")
  # bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  # bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  # coauth<-rbind(coauth,bands_more)
  # 
  # #dim(coauth)
  # new_order<-match(order_Peaks,coauth$name)
  # coauth<-coauth[new_order,]
  # 
  # # grouped
  # #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  # #                                    header = T,sep="\t")
  # Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  # new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  # Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  # 
  # coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  # coauth$grp_mannul<-c(#"DU",
  #                      rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
  #                      rep("Lipids",2),"Lipids; Carbohydrates")
  # coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
  #                                                      "Proteins","Starch","TAG", "Memanbrane Lipids",
  #                                                      "Lipids","Lipids; Carbohydrates"))
  # #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  # #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,
  # #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  # 
  # #coauth$grp_mannul<-coauth$Group_new
  # #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  # #coauth<-coauth[order(coauth$grp_mannul),]
  # 
  # # keep only this people in edges
  # connect <- connect %>%
  #   filter(from %in% coauth$name) %>%
  #   filter(to %in% coauth$name)
  # 
  # # Create a graph object with igraph
  # mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  # 
  # # prepare a vector of n color in the viridis scale
  # mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  # #mycolor <- sample(mycolor, length(mycolor))
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  # #mycolor<-c("#440154ff", "#5cc863ff")
  # #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  # #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  # 
  # #----------------------
  # #no group and assignment
  # p21<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0.4,0), "null"),
  #     panel.spacing=unit(c(0,0,3.4,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  # ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  # 
  # #----------------------
  # #add assignment
  # p22<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
  #        p22,height = 3,width = 3)
  # 
  # #----------------------
  # #add group
  # p23<-ggraph(mygraph2, layout="linear") + 
  #   #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  #   geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
  #   geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
  #   #geom_node_point(aes(size=n+1), alpha=0.5) +
  #   scale_size_continuous(range=c(0.5,3.5)) +
  #   scale_color_manual(values=mycolor) +
  #   #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  #   geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
  #                  colour=mycolor[unclass=coauth$grp_mannul],
  #                  angle=65, 
  #                  hjust=1.1, 
  #                  nudge_y = -1.1, 
  #                  size=1.5) +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0), "null"),
  #     panel.spacing=unit(c(0,0,0,0), "null")
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  # ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
  #        p23,height = 4,width = 2)
  # 
  # 
}
#----------finished-------------



