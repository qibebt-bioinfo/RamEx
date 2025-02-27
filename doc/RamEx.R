## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(RamEx)
library(magrittr)
data(RamEx_data)
options(mc.cores = 2)

## -----------------------------------------------------------------------------
RamEx_data %<>%  Preprocessing.Smooth.Sg %>% Preprocessing.Baseline.Polyfit %>% Preprocessing.Normalize(.,'ch') 
mean.spec(RamEx_data@datasets$normalized.data, RamEx_data@meta.data$group)  

## -----------------------------------------------------------------------------
qc_icod <- Qualitycontrol.ICOD(RamEx_data@datasets$normalized.data,var_tol = 0.5)
data_cleaned <- RamEx_data[qc_icod$quality,] 
mean.spec(data_cleaned@datasets$normalized.data, data_cleaned@meta.data$group,0.3)
qc_mcd <- Qualitycontrol.Mcd(RamEx_data@datasets$normalized.data) 
qc_t2 <- Qualitycontrol.T2(RamEx_data@datasets$normalized.data) 
qc_dis <- Qualitycontrol.Dis(RamEx_data@datasets$normalized.data) 
qc_snr <- Qualitycontrol.Snr(RamEx_data@datasets$normalized.data, 'easy')
## -----------------------------------------------------------------------------
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
# calculate CDR
CDR <- data.frame(data_cleaned@meta.data, 
                  data_cleaned@interested.bands$`2000~2250`/(data_cleaned@interested.bands$`2000~2250` + data_cleaned@interested.bands$`2750~3050`))
## -----------------------------------------------------------------------------
data.reduction <- Feature.Reduction.Pca(data_cleaned, draw=T, save = F) %>% Feature.Reduction.Pcoa(., draw=T, save = F) %>% Feature.Reduction.Tsne(., draw=T, save = F) %>% Feature.Reduction.Umap(., draw=T, save=F) 
## -----------------------------------------------------------------------------
# Markers analysis

ROC_markers <- Raman.Markers.Roc(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],data_cleaned@meta.data$group, paired  = TRUE, threshold = 0.98) 
cor_markers <- Raman.Markers.Correlations(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned@meta.data$group), min.cor = 0.98) 
RBCS.markers <- Raman.Markers.Rbcs(data_cleaned, threshold = 0.003, draw = F) 

# IRCA
#-Global IRCA. This module maybe consumed for a longer period of time due to the image drawing

IRCA.interests <- Intraramanome.Analysis.Irca.Global(data_cleaned)

#-Local IRCA

bands_ann <- data.frame(rbind(cbind(c(742,850,872,971,997,1098,1293,1328,1426,1576),'Nucleic acid'),
                              cbind(c(824,883,1005,1033,1051,1237,1559,1651),'Protein'),
                              cbind(c(1076,1119,1370,2834,2866,2912),'Lipids')))
colnames(bands_ann) <- c('Wave_num', 'Group')
Intraramanome.Analysis.Irca.Local(data_cleaned, bands_ann = bands_ann)

#- 2D-COS

data_cos <- Intraramanome.Analysis.2Dcos(data_cleaned) 

# Phenotype analysis

clusters_louvain <- Phenotype.Analysis.Louvaincluster(object = data_cleaned, resolutions = c(0.8)) 
clusters_kmneans <- Phenotype.Analysis.Kmeans(data_cleaned,5)
clusters_hca <- Phenotype.Analysis.Hca(data_cleaned)


# Classifications
#-PC-LDA
#-SVM
#-Random Forest

model.gmm <- Classification.Gmm(data_cleaned)
model.lda <- Classification.Lda(data_cleaned)
model.rf <- Classification.Rf(data_cleaned)
model.svm <- Classification.Svm(data_cleaned)

# Quantifications

quan_pls <- Quantification.Pls(data_cleaned) 
quan_mlr <- Quantification.Mlr(data_cleaned) 
quan_glm <- Quantification.Glm(data_cleaned) 


# Spectral decomposition

decom_mcr <- Spectral.Decomposition.Mcrals(data_cleaned,2)
decom_ica <- Spectral.Decomposition.Ica(data_cleaned, 2) 
decom_nmf <- Spectral.Decomposition.Nmf(data_cleaned) 


