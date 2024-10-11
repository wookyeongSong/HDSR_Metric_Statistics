# HarvardDataReviewAAL

Welcome to the repository for the Harvard Data Review section, where we conduct a comprehensive analysis using Autism Brain Imaging Data Exchange (ABIDE) dataset, which includes 871 subjects (727 males and 144 females) with ages ranging from 6.47 to 58 years. This dataset was fetched from the Nilearn Python package (https://nilearn.github.io/dev/modules/generated/nilearn.datasets.fetch_abide_pcp.html).

Key details:

Atlas: I selected the AAL Atlas due to its extensive use in Neuroscience and its success in our previous projects. This atlas parcellates the whole brain into 116 ROIs.

Files: 
- aal_labels.csv: Contains the names of the 116 ROIs.
- abide_phenotypic.csv: Includes various covariates such as age, sex, disease status, and site ID. For detailed documentation of these covariates, refer to the link above.
- aal folder: Contains the BOLD time series matrices (196 x 116) for all subjects, with each CSV file named according to the corresponding subject ID.

There are 20 sites in this dataset. To avoid site effects (e.g., different machines), we could focus on the 'NYU' site, which has the largest number of subjects (136 males and 36 females).
