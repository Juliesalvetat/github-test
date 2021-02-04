# Give name to the columns of every csv extracted from Solomon Coder 
# in choosen directory path 

# input : directory set by user 
# output : Processed_GOPR0XXX.csv

# Use ConfigFinalColumnNames to name column with new Solomon configuration 

######################################## Manual user intervention ##########################
# Set directory path 
#setwd("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\FAROFA_2")
setwd("C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\FAROFA_1")
######################################## Manual user intervention ##########################

#install.packages("plyr")
#install.packages("xlsx")
#install.packages("rio")
#install.packages("clipr")
#install.packages("csvy")
#install.packages("feather")
#install.packages("fst")
#install.packages("hexView")
#install.packages("readODS")
#install.packages("rmatio")
#install.packages("xml2")
#install.packages("DescTools")

library("plyr") # for the function rbind.fill
#library("xlsx")
#library("rio")
#library("dplyr")
#library("DescTools")

# Read ConfigFinalColumnNames.csv that contains the column names 
ConfigFinalColumnNames_path<-"C:\\Users\\jsalveta\\Desktop\\processing video\\ConfigFinalColumnNames.csv"
ConfigFinalColumnNames <- read.csv(ConfigFinalColumnNames_path,header = TRUE)

# Find csv extracted from Solomon Coder
# Find all csv that start with letter G in folders and subfolders of directory path 
#csv_from_solomon_path <- list.files(recursive = TRUE, full.name = TRUE)

csv_from_solomon_path <- list.files(pattern = "^G(.*)csv$", recursive = TRUE, full.name = TRUE)

#for (i in 121:126) {
for (i in 1:length(csv_from_solomon_path)) {
  
  # Read SOlomon csv 2 column number and Time
  table_from_solomon <- read.csv(csv_from_solomon_path[i],header = TRUE, sep = ";", quote = "\"", 
                                 dec = ",", fill = FALSE, comment.char="")
  # Define number from Solomon csv as character string 
  number = as.character(table_from_solomon$number)
  # Split number using ; as delimiter into a list
  split_number_list <- strsplit(number, ";")
  # Fill list with NA using length of bigger element in list to omake a list of equal size elements
  split_number_list<-lapply(split_number_list, `length<-`, max(lengths(split_number_list)))
  # Transform that list into a data.frame
  new_table<-rbind.fill(lapply(split_number_list, function(X) data.frame(t(X))))
  
  colnames(new_table) <- colnames(ConfigFinalColumnNames [1:length(new_table)])
  # colnames(new_table) <- ConfigFinalColumnNames_vect[1:length(new_table)] # be carefull to have the same 
  
  # Replace strings 
  new_table[new_table==""]<-NA
  
  new_table <- dplyr::bind_rows(ConfigFinalColumnNames, new_table)
  # Add Time column at the beguining of new_table
  Time <- data.frame(as.character(table_from_solomon$Time))
  colnames(Time)<-"Time"
  new_table <- cbind(Time,new_table)
  
  # Add video informations
  # If your in FAROFA_X folder
  saving_path <- dirname(csv_from_solomon_path[i])
  split_path <- unlist(strsplit(saving_path, "/"))
  
  Cruise <- data.frame(rep(basename(getwd()), nrow(new_table)))
  colnames(Cruise) <- "Cruise"
  VideoType <- data.frame(rep(split_path[3], nrow(new_table)))
  colnames(VideoType) <- "VideoType"
  Date <- data.frame(rep(split_path[4], nrow(new_table)))
  colnames(Date) <- "Date" 
  VideoNumber <- data.frame(rep(split_path[5], nrow(new_table)))
  colnames(VideoNumber) <- "VideoNumber"
  name_with_ext <- basename(csv_from_solomon_path[i])
  name_without_ext <- tools::file_path_sans_ext(name_with_ext)
  VideoName <- data.frame(rep(name_without_ext, nrow(new_table)))
  colnames(VideoName) <- "VideoName"
  new_table <- cbind(Cruise,VideoType,Date,VideoNumber,VideoName, new_table)
  
  # Save new_table in the right place
  new_name <- paste("Processed",name_with_ext, sep="_")
  new_name_and_path <- paste(saving_path,new_name, sep="/")
  
  
  # write csv saving GOPR0XXX_processed in :
  # "FAROFA_X/VIDEOS/VideoType/YYYY-MM-DD/VT_XX/Processed_GOPR0XXX.csv
  write.csv(new_table,new_name_and_path,row.names = FALSE)
  
} # loop end 
