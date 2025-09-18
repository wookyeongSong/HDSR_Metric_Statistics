# HDSR_Metric_Statistics

This repository contains the implementation for the paper "Non-Euclidean Data Analysis With Metric Statistics" which were run using R. This repository contains a comprehensive analysis using Autism Brain Imaging Data Exchange (ABIDE) dataset, fetched from the Nilearn Python package (https://nilearn.github.io/dev/modules/generated/nilearn.datasets.fetch_abide_pcp.html). 

## Dataset
We selected the AAL Atlas due to its extensive use in Neuroscience and its success in our previous projects. This atlas parcellates the whole brain into 116 ROIs. Dataset can be found in `applications/abide_aal` folder.

Files: 
- aal_labels.csv: Contains the names of the 116 ROIs.
- abide_phenotypic.csv: Includes various covariates such as age, sex, disease status, and site ID. For detailed documentation of these covariates, refer to the link above.
- aal folder: Contains the BOLD time series matrices (196 x 116) for all subjects, with each CSV file named according to the corresponding subject ID.

## Details on reproducing results, figures, and tables

In `frechet_regression` folder:
* `/Figures1&2.R`: Creates Figures 1 and 2 demonstrating Fréchet regression on a toy example with weighted networks.

In `variation_analysis` folder:
* `/variance_mainfunctions.R`: Contains all main functions necessary for the simulations and real data analysis regarding variation_analysis in Section 3.
* `/Figure_3.R`: Produces results corresponding to Figure 3 in Section 3.

In `distance_profile` folder:
* `/utlis.R`: Contains all main functions necessary for the distance profile simulations in Section 4.
* `/Figure_4.R`-`/Figure_6.R`: Produce Figures 4, 5, and 6 in Section 4, respectively.

In `applications` folder:
* `/abide_aal` folder: Contains ABIDE dataset.
* `/gnr.R`: Contains the main Fréchet regression function for network data. This implements global network regression using graph Laplacians with Euclidean predictors, supporting both Frobenius and power metrics.
* `/NetworkRegressionABIDE.R`: Implements Fréchet regression analysis for functional connectivity networks using the ABIDE dataset. Compares network properties between autism spectrum disorder (ASD) and control groups across different ages. Produces Figures 7, 8, and 9 showing network visualizations and metrics.
* `/variance_autism.R`: Performs variance analysis of autism brain imaging data for Section 5. Analyzes the variability in functional connectivity networks between ASD and control groups, producing Figure 10.
* `/data_process.R`: Preprocess the data required for two-sample testing in Section 5.
* `/utlis.R`: Contains all main functions necessary for the two sample testing in Section 5.
* `/test.R`: Obtain the p-values of the two sample testing in Section 5.
* `/Figure_10.R`: Produce Figure 10 in Section 5.
* `/Table_1.R`: Produce Table 1 in Section 5.
