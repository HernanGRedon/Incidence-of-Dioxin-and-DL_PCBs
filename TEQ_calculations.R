############################################################################################################
# TEQ calculation
library(readxl)
library(dplyr)
##################
#Load TEF VALUES
TEFvalue <- readxl::read_xlsx("WHO-TEF-2022.xlsx")
TEFvalue <- TEFvalue[,c("Congener","2022 WHO-TEF")]

#Load 
Exposure <- readxl::read_xlsx("ferg2-CTTF-dioxin_BMdata_ALLcompounds_F.xlsx", sheet = "Exposure")


unique_combinations <- Exposure %>%
  distinct(REF_LOCATION, REF_YEAR_START) 


TEQcomp2 <- c()
for( i in 1:nrow(TEFvalue)) {
  
  Compound <- unlist(TEFvalue[i,1])
  TEF <- as.vector(unlist(TEFvalue[i,2]))
  
  TEQcomp <- Exposure %>% filter(SUBSTANCE == Compound)
  
  TEQcomp$VALUE_MEAN_TEQ <-  TEQcomp$VALUE_MEAN * TEF
  TEQcomp$VALUE_MEDIAN_TEQ <-  TEQcomp$VALUE_MEDIAN * TEF
  TEQcomp$VALUE_P000_TEQ <-  TEQcomp$VALUE_P000 * TEF
  TEQcomp$VALUE_P100_TEQ <-  TEQcomp$VALUE_P100 * TEF
  
  TEQcomp2 <- rbind.data.frame(TEQcomp2,TEQcomp)
  
}


TEQcomp2 <- as.data.frame(TEQcomp2)

TEQcomp2$VALUE_UNIT <- "pg TEQ/g fat"
TEQcomp2$VALUE_MEAN <- TEQcomp2$VALUE_MEAN_TEQ
TEQcomp2$VALUE_MEDIAN <- TEQcomp2$VALUE_MEDIAN_TEQ
TEQcomp2$VALUE_P000 <- TEQcomp2$VALUE_P000_TEQ
TEQcomp2$VALUE_P100 <- TEQcomp2$VALUE_P100_TEQ

TEQcomp2$VALUE_MEAN_TEQ <- NULL
TEQcomp2$VALUE_MEDIAN_TEQ<- NULL
TEQcomp2$VALUE_P000_TEQ<- NULL
TEQcomp2$VALUE_P100_TEQ<- NULL


writexl::write_xlsx(TEQcomp2,"ferg2-CTTF-dioxin-TEQ-06_03_2024.xlsx")
#######################################################################################