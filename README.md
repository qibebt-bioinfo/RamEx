<p align="center">
  <a href="https://github.com/qibebt-bioinfo/RamEx">
    <img src="https://ramex-1257170263.cos.ap-guangzhou.myqcloud.com/RamEx-logo.jpg" alt="RamEx logo"  width="300">
  </a>
</p>

## *RamEx*: An R package for high-throughput microbial ramanome analyses with accurate quality assessment

## Key features

![Overview of RamEx](https://ramex-1257170263.cos.ap-guangzhou.myqcloud.com/RamEx-pipeline.jpg)
- **Reliability** achieved via stringent statistical control
- **Robustness** achieved via flexible modelling of the data and automatic parameter selection
- **Reproducibility** promoted by thorough recording of all analysis steps
- **Ease of use**: high degree of automation, an analysis can be set up in several mouse clicks, no bioinformatics expertise required
- **Powerful tuning options** to enable unconventional experiments
- **Scalability and speed**: up to 100 runs processed per minutes

**Download**: https://github.com/qibebt-bioinfo/RamEx


### Table of Contents
**[Installation](#installation)**<br>
**[Getting started](#getting-started)**<br>
**[Raw data formats](#raw-data-formats)**<br>
**[Output](#output)**<br>
**[Changing default settings](#changing-default-settings)**<br>
**[Visualisation](#visualisation)**<br>
**[Frequently asked questions (FAQ)](#frequently-asked-questions)**<br>
**[Contact](#contact)**<br>

### Installation

RamEx will be installed from GitHub:.
```{r}
library('devtools')
install_github("qibebt-bioinfo/RamEx")
```

### Getting Started
#### Data Loading
Raman spectra are respectively tracked in single txt files, and their meta info is recorded in the file name.
Here we assume there's only one factor of the dataset, which means RamEx do not contain multiple-factor analysis. If you have multiple factors but they are independent of each other, these factors will be treated as one factor.
```{r}
library(RamEx)
library(magrittr)
data(RamEx_data)
```
#### Pretreatment
Spectral pretreatment will make the spectrum clearer, containing smoothing, baseline removal, normalization and truncation.
Mean spectra will display their effects.
Here the results of each step will be kept in the Ramanome for better debugging, and 'draw.mean' exhibit the final dataset.
```{r}
RamEx_data %<>%  Preprocessing.Smooth.Sg %>% Preprocessing.Baseline.Polyfit %>% Preprocessing.Normalize(.,'ch') 
mean.spec(RamEx_data@datasets$normalized.data, RamEx_data@meta.data$group)  
```
#### Quality control
```{r}
qc_icod <- Qualitycontrol.ICOD(RamEx_data,var_tol = 0.5)
data_cleaned <- RamEx_data[qc_icod$quality,] 
mean.spec(data_cleaned@datasets$normalized.data, data_cleaned@meta.data$group,0.3)
qc_mcd <- Qualitycontrol.Mcd(RamEx_data) 
qc_t2 <- Qualitycontrol.T2(RamEx_data) 
qc_dis <- Qualitycontrol.Dis(RamEx_data) 
qc_snr <- Qualitycontrol.Snr(RamEx_data, 'easy')
```

#### Interested Bands
Get single-cell intensitiy or intensity accumulationy within a wavenumber range, pls give a list containing multiple bands or band ranges. These feature selection results will be saved as 'interested.bands' in the given Ramanome object. Further, you can add some equations by yourself. 
```{r}
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
# calculate CDR
CDR <- data.frame(data_cleaned@meta.data, 
                  data_cleaned@interested.bands$`2000~2250`/(data_cleaned@interested.bands$`2000~2250` + data_cleaned@interested.bands$`2750~3050`))
```

#### Reduction
Nonlinear methods, such as UMAP and t-SNE. Linear methods like PCA, pCoA. The reduced sample matrix will be contained in the Ramanome onject as 'reductions'. Attention: RamEx uses PCA to reduce the dimensions of the high-dimensional spectrum, since UMAP and t-SNE are highly complex algorithms.
```{r}
data.reduction <- Feature.Reduction.Pca(data_cleaned, draw=T, save = F) %>% Feature.Reduction.Pcoa(., draw=T, save = F) %>% Feature.Reduction.Tsne(., draw=T, save = F) %>% Feature.Reduction.Umap(., draw=T, save=F) 
``` 

#### Markers analysis
```{r}
ROC_markers <- Raman.Markers.Roc(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],data_cleaned@meta.data$group, paired  = TRUE, threshold = 0.8) 
cor_markers <- Raman.Markers.Correlations(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned@meta.data$group), min.cor = 0.8) 
RBCS.markers <- Raman.Markers.Rbcs(data_cleaned, threshold = 0.003, draw = F) 
```
#### IRCA
-Global IRCA. This module maybe consumed for a longer period of time due to the image drawing
```{r}
IRCA.interests <- Intraramanome.Analysis.Irca.Global(data_cleaned)
```
-Local IRCA
```{r}
bands_ann <- data.frame(rbind(cbind(c(742,850,872,971,997,1098,1293,1328,1426,1576),'Nucleic acid'),
                              cbind(c(824,883,1005,1033,1051,1237,1559,1651),'Protein'),
                              cbind(c(1076,1119,1370,2834,2866,2912),'Lipids')))
colnames(bands_ann) <- c('Wave_num', 'Group')
Intraramanome.Analysis.Irca.Local(data_cleaned, bands_ann = bands_ann)
```
- 2D-COS
```{r}
data_cos <- Intraramanome.Analysis.2Dcos(data_cleaned) 
```
#### Phenotype analysis
```{r}
clusters_louvain <- Phenotype.Analysis.Louvaincluster(object = data_cleaned, resolutions = c(0.8)) 
clusters_kmneans <- Phenotype.Analysis.Kmeans(data_cleaned,5)
clusters_hca <- Phenotype.Analysis.Hca(data_cleaned)
```

#### Classifications
-PC-LDA
-SVM
-Random Forest
```{r}
model.gmm <- Classification.Gmm(data_cleaned)
model.lda <- Classification.Lda(data_cleaned)
model.rf <- Classification.Rf(data_cleaned)
model.svm <- Classification.Svm(data_cleaned)
```
#### Quantifications
```{r}
quan_pls <- Quantification.Pls(data_cleaned) 
quan_mlr <- Quantification.Mlr(data_cleaned) 
quan_glm <- Quantification.Glm(data_cleaned) 
```

#### Spectral decomposition
```{r}
decom_mcr <- Spectral.Decomposition.Mcrals(data_cleaned,2)
decom_ica <- Spectral.Decomposition.Ica(data_cleaned, 2) 
decom_nmf <- Spectral.Decomposition.Nmf(data_cleaned) 
```
### Raw data formats

It accommodates data from mainstream instrument manufactures such as Horiba, Renishaw, Thermo Fisher Scientific, WITec, and Bruker. This module efficiently manages single-point data collection, where each spectrum is stored in a separate txt file, as well as mapping data enriched with coordinate information. 

### Output

The **Output** pane allows to specify where the output should be saved. 

### Changing default settings
RamEx can be successfully used to process almost any experiment with default settings. In general, it is recommended to only change settings when specifically advised to help information.


### Visualisation
RamEx also offers an online version. Please visit (http://ramex.single-cell.cn).


<!--### functions

**Import** Import the xx data
* **R Data** xxx
**Quality Control** Import the xx data
* **Outlier Detection** xxx
**Cell-level analysis** Import the xx data
* **Outlier Detection** xxx
**Singal-level analysis** Import the xx data
* **RBCS** xxx
**Visualization** Import the xx data
* **mean_spectrum** xxx---> 

### Frequently asked questions
**Q: Why RamEx?**  
**A:** Raman spectroscopy, with its fast, label-free, and non-destructive nature, is increasingly popular for capturing vibrational energy levels and metabolic differences in cells, providing qualitative and quantitative insights at single-cell or subcellular resolutions. Leveraging the extensive information provided by the complex and high-dimensional nature of Ramanome, we developed RamEx, an R package designed to adeptly manage extensive Raman datasets generated by a wide range of devices and instruments. It features: 1) a dynamic outlier detection algorithm that operates without prior knowledge or fixed criteria; 2) optimized clustering and marker identification algorithms tailored to the unique properties of high dimensional, colinear and nonlinear Raman spectra; 3 ) a unified computational framework with tools and pipelines for key Raman tasks such as cell type/species identification, clusteringphenotypic analysis, and antibiotic resistance detectionmolecular composition analysis; 4) enhanced processing of large-scale datasets through C++ optimization and GPU computing; 5) a standardized Raman dataset format with integrated metadata and evaluation metrics; and 6) a graphical user-interface (GUI) for intuitive data visualization and interaction.

**Q: How to install OpenCL?**  
**A:** RamEx relies on OpenCL for GPU acceleration. If you encounter issues with OpenCL installation or configuration, please refer to the following additional resources:
- Linux:
  - Ensure you have the latest drivers for your GPU.
  - For more information on OpenCL on Linux, refer to the [POCL documentation](https://github.com/pocl/pocl).
- Windows:
  - Intel and AMD provide detailed installation guides for their OpenCL SDKs. Follow the instructions specific to your hardware:
    - [Intel OpenCL SDK](https://www.intel.com/content/www/us/en/developer/tools/opencl-sdk.html)
    - [AMD OpenCL SDK](https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases)
- macOS:
  - macOS users can leverage Homebrew for a seamless installation experience also refer to the [POCL documentation](https://github.com/pocl/pocl).
If you still face issues, please consider reaching out to the respective SDK support forums or the [OpenCL community](https://www.khronos.org/opencl/).

<!--(it's recommended to use the latest version - RamEx 2.1)  -->

<!--<img src="xxx"></br> 

<!--Please cite   
**RamEx : xxxx 

[xxxx, 2024](https://github.com/qibebt-bioinfo/RamEx)

When using RamEx for xxxx  ---> 

### Key papers 
**IRCA**   
He Y., Huang S., Zhang P., Ji Y., Xu J., 2021. [Intra-Ramanome Correlation Analysis Unveils Metabolite Conversion Network from an Isogenic Population of Cells](https://doi.org/10.1128/mbio.01470-21). *mBio* 

**RBCS**  
Teng L., Wang X.,  Wang X.,  Gou H.,  Ren L., & Wang T., 2016. [Label-free, rapid and quantitative phenotyping of stress response in e. coli via ramanome](https://www.nature.com/articles/srep34359.pdf). *Scientific Reports* 

<!--**Other key papers**  
- IRCA:  
He Y,,Huang S,,,Zhang P,,Ji Y, Xu J,, 2021. Intra-Ramanome Correlation Analysis Unveils Metabolite Conversion Network from an Isogenic Population of Cells. mBio 12:10.1128/mbio.01470-21.
- xxxx:   
[x'x'x, 2021](https://github.com/qibebt-bioinfo/RamEx)

**R package** some useful functions:https://github.com/qibebt-bioinfo/RamEx 
**Visualisation** of ramanome-->




### Contact
Please post any questions, feedback, comments or suggestions on the [GitHub Discussion board](https://github.com/qibebt-bioinfo) (alternatively can create a GitHub issue) or email SCC.
