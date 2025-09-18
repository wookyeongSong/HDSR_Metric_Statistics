# HDSR_Metric_Statistics

This repository contains the implementation for the paper ``Non-Euclidean Data Analysis With Metric Statistics" which were run using R. This repository contains a comprehensive analysis using Autism Brain Imaging Data Exchange (ABIDE) dataset. This dataset was fetched from the Nilearn Python package (https://nilearn.github.io/dev/modules/generated/nilearn.datasets.fetch_abide_pcp.html).

## Details on reproducing results, figures, and tables

In `frechet_regression` folder:
* To be updated
* To be updated

In `variation_analysis` folder:
* `/HDSR_Variance_mainfunctions.R`: Contains all main functions necessary for the simulations and real data analysis regarding variation_analysis.
* `/HDSR_Variance_Sphere.R`: Produces results corresponding to Figure 3 in Section 3.
* `/HDSR_Variance_Autism.R`: Produces results corresponding to Table 1 and Figure 10 in Section 5.

In `variation_analysis` folder:
* `/HDSR_Variance_mainfunctions.R`: Contains all main functions necessary for the simulations and real data analysis regarding variation_analysis.
* `/HDSR_Variance_Sphere.R`: Produces results corresponding to Figure 3 in Section 3.
* `/HDSR_Variance_Autism.R`: Produces results corresponding to Table 1 and Figure 10 in Section 5.
  
Atlas: I selected the AAL Atlas due to its extensive use in Neuroscience and its success in our previous projects. This atlas parcellates the whole brain into 116 ROIs.

Files: 
- aal_labels.csv: Contains the names of the 116 ROIs.
- abide_phenotypic.csv: Includes various covariates such as age, sex, disease status, and site ID. For detailed documentation of these covariates, refer to the link above.
- aal folder: Contains the BOLD time series matrices (196 x 116) for all subjects, with each CSV file named according to the corresponding subject ID.

There are 20 sites in this dataset. To avoid site effects (e.g., different machines), we could focus on the 'NYU' site, which has the largest number of subjects (136 males and 36 females).
