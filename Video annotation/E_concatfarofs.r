library(lubridate)
library(R.matlab)

paht2data <- "C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\Big_table_FAROFA_1_coordinates.csv"
data_video <- read.csv(paht2data)
data_video <- data_video[which(data_video$Situation=="in"),]
big_table <- data_video

# convert date field to R date format
data_video$Date = as.character(data_video$Date)
x <- paste(data_video$Date, data_video$Time_UTC)
date_video_formatR <- strptime(x, "%Y-%m-%d %HH %MM %OS",tz='UTC')
x[which(is.na(date_video_formatR))]
#data_video$Time_UTC[which(is.na(date_video_formatR))]
#unique(data_video$VideoNumber[which(is.na(date_video_formatR))])
#tmp<-data_video[which(is.na(date_video_formatR)),]
#tmp<-data_video$Time_UTC[which(is.na(date_video_formatR)),]
#data_video$VideoNumber[which(is.na(date_video_formatR))]

# reading EK80 echointegration file  
data_EK80_FM <- readMat('C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie/FAROFA1_FM.mat')

# Matlab dates conversion into R date format
date_EK80_formatR_FM = as.POSIXct((data_EK80_FM$Time-719529)*86400,origin="1970-01-01",tz='UTC')
# For each GoPro observation, looking for the closest lon/lat from EK80 file

data_EK80 <- readMat('C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie/FAROFA1.mat')
date_EK80_formatR = as.POSIXct((data_EK80$Time-719529)*86400,origin="1970-01-01",tz='UTC')
data_video<-data_video[which(is.na(data_video$Latitude)),]

for (i in 1:length(date_video_formatR)) {
  i
  ind = min(which(abs(as.numeric(difftime(date_video_formatR[i], date_EK80_formatR, units = "secs")))<seconds(1)))
  data_video$Longitude[i] = data_EK80$Longitude[ind];
  data_video$Latitude[i] = data_EK80$Latitude[ind];
  data_video$Bottom[i] = data_EK80$Bottom[ind];
}

paht2save <- "C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\Big_table_FAROFA_1_coordinates_test1.csv"
write.csv(data_video,paht2save,row.names = FALSE)



library(plyr)
setwd("C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\25-09-2020")

new_name_and_path<-"C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\25-09-2020\\final_25-09-2020.csv"
# list all fils with name starting with "Processed"
csv_processed_path <- list.files(pattern = "*coordinates.csv", recursive = TRUE, full.name = TRUE)

# read all processed csv 
datalist= lapply(csv_processed_path, function (x) read.csv(file=x, sep = ",", header=TRUE,colClasses=c("character")))
big_table <- do.call(rbind.fill, lapply(datalist, function(x) x[match(names(as.data.frame(datalist[1])), names(x))]))
tmp <- big_table[which(big_table$Situation=="in"),]
tmp2 <- tmp[which(is.na(big_table$Latitude)),]
a<-cbind(tmp[c(1:3,74:76,4:6,73,8)])

b<-tmp[c(9:72)]
d<-cbind(b[c(1:8,39,41,47,60:63)])
oldcolnames <- colnames(d)
newcolnames <- c("AlgaePelagicGreen","AlgaePelagicBrown","AlgaePelagicRed","Gelatinous","Salps","Unknown","FishUn","SharkUn","AlgaePelagicUn","GelatinousNum","RayUn","FishPelagicSmall","CarangidaeSmall","FishReef","FishBottomSmall")
colnames(d) <- newcolnames
d<-  d[ , order(names(d))]
b<-cbind(b[c(9:38,40,42:46,48:59,64)])
b<-  b[ , order(names(b))]
c<-cbind(a,b,d)
write.csv(c,new_name_and_path,row.names = FALSE)
plot(c$Longitude,c$Latitude)

unique(c$VideoNumber[which(is.na(c$Latitude))])
