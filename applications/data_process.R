source("utlis.R")
library(ODP)
library(energy)
data_all_chorts <- read.csv("abide_aal/abide_phenotypic.csv")
cohort_names<-unique(data_all_chorts$SITE_ID)
for (cohort_Int in cohort_names) {
  data_all <- read.csv("abide_aal/abide_phenotypic.csv")
  data_Int<-data_all[data_all$SITE_ID==cohort_Int,]
  
  
  Data <- list() 
  count <- 1      
  for (id in data_Int$subject) {
    # Use tryCatch to handle warnings
    # Read the file for the given id
    temp <- read.csv(paste0("abide_aal/aal/", id, ".csv"),header = FALSE)
    
    # Check for constant columns (those with zero standard deviation)
    if (!any(apply(temp, 2, sd) == 0)) {
      Data[[count]] <- list(
        id = id,
        mat = cor(temp),
        sex = data_Int[data_Int$subject == id, ]$SEX,
        age = data_Int[data_Int$subject == id, ]$AGE_AT_SCAN,
        autism=data_Int[data_Int$subject == id, ]$DX_GROUP,
        DSM=data_Int[data_Int$subject == id, ]$DSM_IV_TR
      )
      count <- count + 1  # Increment counter
    }
  }
  save(Data,file = paste0(cohort_Int,".RData"))
  
}



