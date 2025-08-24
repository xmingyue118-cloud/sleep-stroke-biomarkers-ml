# sleep-stroke-biomarkers-ml
Code and data for the paper [Exploring common circulating diagnostic biomarkers for sleep disorders and stroke based on machine learning]

## Data
- **Source**: The raw data are publicly available from the GEO database under accession numbers `GSE208668`,`GSE16561`,`GSE22255`and `GSE98566`

## Instructions
1.  **Install R** (version ≥ 4.1.0) and the required packages by running `install.packages(c('tidyverse', 'limma', 'randomForest', 'xgboost', 'pROC'))`.
2.  **Run the scripts in order**:
    - `0_Dataprogress.R`
    - `1_DEGs_WGCNA_GSEA.R`
    - `2_Ai_ROC.R` (Includes feature selection & machine learning model building)
## Notes
- The code has been tested on [Windows 11] with R version [e.g., 4.4.1].
- For questions, please open an Issue here.
