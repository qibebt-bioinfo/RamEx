## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(RamEx)
library(dplyr)
data(RamEx_data)
options(mc.cores = 2)

## -----------------------------------------------------------------------------
RamEx_data <- RamEx_data %>% Preprocessing.Smooth.Sg %>% Preprocessing.Baseline.Polyfit %>% Preprocessing.Normalize(.,'ch') 
mean.spec(RamEx_data$normalized.data, RamEx_data$group)  

## -----------------------------------------------------------------------------
qc_icod <- Qualitycontrol.ICOD(RamEx_data)
data_cleaned <- RamEx_data[qc_icod$quality,] 
mean.spec(data_cleaned$normalized.data, data_cleaned$group,0.3)
qc_mcd <- Qualitycontrol.Mcd(RamEx_data) 
qc_t2 <- Qualitycontrol.T2(RamEx_data) 
qc_dis <- Qualitycontrol.Dis(RamEx_data) 
qc_snr <- Qualitycontrol.Snr(RamEx_data, 'easy')
## -----------------------------------------------------------------------------
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
# calculate CDR
CDR <- data.frame(data_cleaned@meta.data, 
                  data_cleaned$`2000~2250`/(data_cleaned$`2000~2250` + data_cleaned$`2750~3050`))
## -----------------------------------------------------------------------------
data.reduction <- Feature.Reduction.Pca(data_cleaned) %>% Feature.Reduction.Pcoa %>% Feature.Reduction.Tsne %>% Feature.Reduction.Umap
## -----------------------------------------------------------------------------
# Markers analysis

ROC_markers <- Raman.Markers.Roc(Preprocessing.Cutoff(data_cleaned, 1400, 1500), paired  = TRUE, threshold = 0.6) 
cor_markers <- Raman.Markers.Correlations(Preprocessing.Cutoff(data_cleaned, 1400, 1500), min.cor = 0.6) 
RBCS.markers <- Raman.Markers.Rbcs(data_cleaned, threshold = 0.003, show = F) 

# IRCA
#-Global IRCA. This module maybe consumed for a longer period of time due to the image drawing

IRCA.interests <- Intraramanome.Analysis.Irca.Global(data_cleaned)

#-Local IRCA

bands_ann <- data.frame(
  Wave_num = c(742,850,872,971,997,1098,1293,1328,1426,1576,
               824,883,1005,1033,1051,1237,1559,1651,
               1076,1119,1370,2834,2866,2912),
  Group = rep(c('Nucleic acid', 'Protein', 'Lipids'), times = c(10, 8, 6)) )
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
pred_new <- predict_classification(model.lda, data_cleaned)

# Quantifications

quan_pls <- Quantification.Pls(data_cleaned) 
quan_mlr <- Quantification.Mlr(data_cleaned) 
quan_glm <- Quantification.Glm(data_cleaned) 
pred_new <- predict_quantification(quan_pls, data_cleaned)

# Spectral decomposition

decom_mcr <- Spectral.Decomposition.Mcrals(data_cleaned, n_comp = 2)
decom_ica <- Spectral.Decomposition.Ica(data_cleaned, n_comp = 2) 
decom_nmf <- Spectral.Decomposition.Nmf(data_cleaned, n_comp = 2) 


