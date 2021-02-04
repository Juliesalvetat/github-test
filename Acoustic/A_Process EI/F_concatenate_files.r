setwd("C:\\Users\\jsalveta\\Desktop\\Data\\EI FISH NO FISH")
library("plyr") 


csv_path <- list.files(pattern = "Average(.*)csv", recursive = TRUE, full.name = TRUE)
datalist= lapply(csv_path, function (x) read.csv(file=x, sep = ",", header=TRUE))
# bind all table together 
my.dataframe<- do.call(rbind.fill, lapply(datalist, function(x) x[match(names(as.data.frame(datalist[1])), names(x))]))
write.csv(my.dataframe,'F123_Averaged_SaSum.csv',row.names = FALSE)
