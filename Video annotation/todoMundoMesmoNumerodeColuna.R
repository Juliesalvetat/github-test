#install.packages("pkgconfig")
library(dplyr)
library(pkgconfig)
setwd("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\FAROFA_1")

#csv_processed_path <- list.files(pattern = "^Processed_(.*)csv", recursive = TRUE, full.name = TRUE)
csv_processed_path <- list.files(pattern = "^With_Time_Processed_(.*)csv", recursive = TRUE, full.name = TRUE)

ConfigFinalColumnNames_path<-"C:\\Users\\jsalveta\\Desktop\\processing video\\ConfigFinalColumnNames.csv"
ConfigFinalColumnNames <- read.csv(ConfigFinalColumnNames_path,header = TRUE)


for (i in 1:length(csv_processed_path)) {
#for (i in 121:126) {
#for (i in 41:45) {  
  # Read SOlomon csv 2 column number and Time
  processed_table <- read.csv(csv_processed_path[i])
  #names(processed_table)[names(processed_table) == 'trachurus.spp'] <- 'dec.mac'
  #names(processed_table)[names(processed_table) == 'thunnus.sp'] <- 'thunnus.spp'
  #names(processed_table)[names(processed_table) == 'caranx.sp'] <- 'caranx.spp'
  #names(processed_table)[names(processed_table) == 'bothidae.sp'] <- 'bothidae.spp'
  #names(processed_table)[names(processed_table) == 'gobiidae.sp'] <- 'gobiidae.spp'
  #names(processed_table)[names(processed_table) == 'ostraciidae.sp'] <- 'ostraciidae.spp'
  names(processed_table)[names(processed_table) == 'rayUn'] <- 'RayUn'
  names(processed_table)[names(processed_table) == 'GelatinousUn'] <- 'GelatinousNum'
  names(processed_table)[names(processed_table) == 'bothidae_sp'] <- 'bothidae_spp'
  names(processed_table)[names(processed_table) == 'gobiidae_sp'] <- 'gobiidae_spp'
  names(processed_table)[names(processed_table) == 'ostraciidae_sp'] <- 'ostraciidae_spp'
  names(processed_table)[names(processed_table) == 'trachurus_spp'] <- 'dec_mac'
  names(processed_table)[names(processed_table) == 'trachurus_sp'] <- 'dec_mac'
  names(processed_table)[names(processed_table) == 'thunnus_sp'] <- 'thunnus_spp'
  names(processed_table)[names(processed_table) == 'caranx_sp'] <- 'caranx_spp'
  
#  oldcolnames <- colnames(processed_table)
#  newcolnames <-gsub("_",".",oldcolnames)
#  colnames(processed_table) <- newcolnames
  df <- data.frame(matrix(ncol = 6, nrow = 0)) 
  colnames(df) <- names(processed_table[1:6])
 FinalColumnNames <- cbind(df,ConfigFinalColumnNames)
   oldcolnames <- colnames(FinalColumnNames)
    newcolnames <-gsub("\\.","_",oldcolnames)
  colnames(FinalColumnNames) <- newcolnames
 processed_table <- dplyr::bind_rows(FinalColumnNames, processed_table)
 # colnames(processed_table)
  
   write.csv(processed_table ,csv_processed_path[i],row.names = FALSE)
}
