## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(RamEx)
library(magrittr)
data(RamEx_data)
data <- RamEx_data
options(mc.cores = 2)

## -----------------------------------------------------------------------------
data_spike <- Preprocesssing.Background.Spike(data, "CPU") 
data_smoothed <- Preprocesssing.Smooth.Sg(data) 
data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
Preprocesssing.Cutoff(data_normalized,550, 1800)
mean.spec(data_normalized@datasets$baseline.data, data@meta.data$group) 

## -----------------------------------------------------------------------------
data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
data_cleaned <- data_normalized[data_cleaned$index_good,] 
mean.spec(data_cleaned@datasets$normalized.data, data_cleaned@meta.data$group)
#qc_icod <- Qualitycontrol.Mcd(data_normalized@datasets$normalized.data)
#qc_t2 <- Qualitycontrol.T2(data_normalized@datasets$normalized.data) 
qc_dis <- Qualitycontrol.Dis(data_normalized@datasets$normalized.data)
hist(qc_dis$dis)
qc_snr <- Qualitycontrol.Snr(data_normalized@datasets$normalized.data) 

## -----------------------------------------------------------------------------
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
# calculate CDR
CDR <- data.frame(data_cleaned@meta.data, 
                  data_cleaned@interested.bands$`2000~2250`/(data_cleaned@interested.bands$`2000~2250` + data_cleaned@interested.bands$`2750~3050`))

## -----------------------------------------------------------------------------
data_cleaned <- Feature.Reduction.Pca(data_cleaned, draw=TRUE, save = FALSE) %>% Feature.Reduction.Pcoa(., draw=TRUE, save = FALSE) %>% Feature.Reduction.Tsne(., draw=TRUE, save = FALSE) %>% Feature.Reduction.Umap(., draw=TRUE, save=FALSE) #

## -----------------------------------------------------------------------------
ROC_markers <- Raman.Markers.Roc(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],data_cleaned@meta.data$group) 
#cor_markers <- Raman.Markers.Correlations(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned@meta.data$group), min.cor = 0.6)
RBCS.markers <- Raman.Markers.Rbcs(data_cleaned, threshold = 0.003, draw = FALSE) 

## -----------------------------------------------------------------------------
IRCA.interests <- Intraramanome.Analysis.Irca.Global(data_cleaned)

## -----------------------------------------------------------------------------
bands_ann <- data.frame(rbind(cbind(c(742,850,872,971,997,1098,1293,1328,1426,1576),'Nucleic acid'),
                              cbind(c(824,883,1005,1033,1051,1237,1559,1651),'Protein'),
                              cbind(c(1076,1119,1370,2834,2866,2912),'Lipids')))
colnames(bands_ann) <- c('Wave_num', 'Group')
Intraramanome.Analysis.Irca.Local(data_cleaned, bands_ann = bands_ann)

## -----------------------------------------------------------------------------
#Intraramanome.Analysis.2Dcos(data_cleaned) 

## -----------------------------------------------------------------------------
#clusters_louvain <- Phenotype.Analysis.Louvaincluster(object = data_cleaned, resolutions = c(0.8)) 
clusters_kmneans <- Phenotype.Analysis.Kmeans(data_cleaned)
clusters_hca <- Phenotype.Analysis.Hca(data_cleaned) 

## -----------------------------------------------------------------------------
Classification.Gmm(data_cleaned) 
Classification.Lda(data_cleaned)
Classification.Rf(data_cleaned)
Classification.Svm(data_cleaned)

## -----------------------------------------------------------------------------
quan_pls <- Quantification.Pls(data_cleaned)
quan_mlr <- Quantification.Mlr(data_cleaned)
quan_glm <- Quantification.Glm(data_cleaned)

## -----------------------------------------------------------------------------
decom_mcr <- Spectral.Decomposition.Mcrals(data_cleaned,2)
decom_ica <- Spectral.Decomposition.Ica(data_cleaned, 2) 
data_nmf <- data_cleaned
data_nmf@datasets$normalized.data %<>% abs
decom_nmf <- Spectral.Decomposition.Nmf(data_nmf)

