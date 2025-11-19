<p align="center">
  <a href="https://github.com/qibebt-bioinfo/RamEx">
    <img src="doc/RamEx-LOGO.jpg" alt="RamEx logo" width="300">
  </a>
</p>

## *RamEx*: An R package for high-throughput microbial ramanome analyses with accurate quality assessment

## Key features

![Overview of RamEx](doc/OverView.jpg)
- **Reliability** achieved via stringent statistical control
- **Robustness** achieved via flexible modelling of the data and automatic parameter selection
- **Reproducibility** promoted by thorough recording of all analysis steps
- **Ease of use**: high degree of automation, an analysis can be set up in several mouse clicks, no bioinformatics expertise required
- **Powerful tuning options** to enable unconventional experiments
- **Scalability and speed**: up to 100 runs processed per minutes

**More Details**: ![Wiki] (https://github.com/qibebt-bioinfo/RamEx/wiki)
**Usage**: [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/qibebt-bioinfo/RamEx)

## Quick Start

### Installation
```
library('devtools')
install_github("qibebt-bioinfo/RamEx")
```
### Data Loading
```{r}
library(RamEx)
library(magrittr)
data(RamEx_data)
```

### Pretreatment
```{r}
RamEx_data %<>%  Preprocessing.Smooth.Sg %>% Preprocessing.Baseline.Polyfit %>% Preprocessing.Normalize(.,'ch') 
mean.spec(RamEx_data@datasets$normalized.data, RamEx_data@meta.data$group)  
```

### Quality control
```{r}
qc_icod <- Qualitycontrol.ICOD(RamEx_data@datasets$normalized.data,var_tol = 0.5)
data_cleaned <- RamEx_data[qc_icod$quality,] 
mean.spec(data_cleaned@datasets$normalized.data, data_cleaned@meta.data$group,0.3)
```

### Interested Bands
```{r}
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
# calculate CDR
CDR <- data.frame(data_cleaned@meta.data,
                  data_cleaned$`2000~2250`/(data_cleaned$`2000~2250` + data_cleaned$`2750~3050`))
```

### Reduction
```{r}
data.reduction <- Feature.Reduction.Umap(data_cleaned, draw=T, save = F)
``` 

### Markers analysis
```{r}
ROC_markers <- Raman.Markers.Roc(data_cleaned$normalized.data[,sample(1:1000, 50)],data_cleaned$group, paired  = TRUE, threshold = 0.8) 
cor_markers <- Raman.Markers.Correlations(data_cleaned$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned$group), min.cor = 0.8) 
```

### IRCA
```{r}
bands_ann <- data.frame(rbind(cbind(c(742,850,872,971,997,1098,1293,1328,1426,1576),'Nucleic acid'),
                              cbind(c(824,883,1005,1033,1051,1237,1559,1651),'Protein'),
                              cbind(c(1076,1119,1370,2834,2866,2912),'Lipids')))
colnames(bands_ann) <- c('Wave_num', 'Group')
Intraramanome.Analysis.Irca.Local(data_cleaned, bands_ann = bands_ann)
```


### Phenotype analysis
```{r}
clusters_louvain <- Phenotype.Analysis.Louvain(object = data_cleaned, resolutions = c(0.8)) 
```

### Classifications
```{r}
model.svm <- Classification.Svm(data_cleaned)
```

### Quantifications
```{r}
quan_pls <- Quantification.Pls(data_cleaned) 
```

### Spectral decomposition
```{r}
decom_mcr <- Spectral.Decomposition.Mcrals(data_cleaned,2)
```

## Raw data formats
It accommodates data from mainstream instrument manufactures such as Horiba, Renishaw, Thermo Fisher Scientific, WITec, and Bruker. This module efficiently manages single-point data collection, where each spectrum is stored in a separate txt file, as well as mapping data enriched with coordinate information. 
Specifically, RamEx automatically interpret and reshape various spectral data layouts, supporting both row‑wise and column‑wise sample arrangements.  
File structure detection rules:
  Type 1: Two columns (wavenumber, intensity)
  Type 2: Mapping matrix — first column is wavenumber, remaining columns are multiple spectra
  Type 3: Coordinate scan — columns 1–n are (x1, ..., xn) metadata annotations, column n+1 is wavenumber, column n+2 is intensity
The above applies to the sample of rows as well.

## Key papers 
**RamEx**   
Zhang Y., Jing G., ..., Xu J., Sun L., 2025. [RamEx: An R package for high-throughput microbial ramanome analyses with accurate quality assessment](https://doi.org/10.1101/2025.03.10.642505). *bioRxiv* 

**IRCA**   
He Y., Huang S., Zhang P., Ji Y., Xu J., 2021. [Intra-Ramanome Correlation Analysis Unveils Metabolite Conversion Network from an Isogenic Population of Cells](https://doi.org/10.1128/mbio.01470-21). *mBio* 

**RBCS**  
Teng L., ...,  Huang W.E., Xu J., 2016. [Label-free, rapid and quantitative phenotyping of stress response in e. coli via ramanome](https://www.nature.com/articles/srep34359.pdf). *Scientific Reports* 



### Contact
Please post any questions, feedback, comments or suggestions on the [GitHub Discussion board](https://github.com/qibebt-bioinfo) (alternatively can create a GitHub issue) or email SCC.
