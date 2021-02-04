library("plyr") 
######################################## Manual user intervention ##########################
ConfigFinalColumnNames_path<-"C:\\Users\\jsalveta\\Desktop\\processing video\\ConfigFinalColumnNames.csv"
ConfigFinalColumnNames <- read.csv(ConfigFinalColumnNames_path,header = TRUE)
# Set directory path 
setwd("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\FAROFA_2\\")
######################################## Manual user intervention ##########################
setwd("C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\FAROFA_2\\")

# list all fils with name starting with "Processed"
csv_processed_path <- list.files(pattern = "^With_Time_Processed_(.*)csv", recursive = TRUE, full.name = TRUE)

# read all processed csv 
datalist= lapply(csv_processed_path, function (x) read.csv(file=x, sep = ",", header=TRUE,colClasses=c("character")))


# bind all table together 
big_table <- do.call(rbind.fill, lapply(datalist, function(x) x[match(names(as.data.frame(datalist[1])), names(x))]))

# Save big_table in the right place
saving_path <- dirname(getwd())
Cruise <- basename(getwd())
new_name <- paste(Cruise,"Big_table.csv", sep="/")
new_name <- paste("Big_table",Cruise, sep="_")
new_name <- paste(new_name,".csv", sep="")
new_name_and_path <- paste(saving_path,new_name, sep="/")

# write csv saving big_table in path set by user at the beguining of the code
write.csv(big_table,new_name_and_path,row.names = FALSE)

