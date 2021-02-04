setwd("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie")
library(lubridate)
library(ggmap)

paht2data <- "C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\24-11-2020\\FAROFA_underwater_video_species_observations.csv"
data_video_all <- read.csv(paht2data,header = TRUE, sep = ";")

#select DV
data_video = data_video_all[which(data_video_all$VideoType=="TowedVideos"),]
# convert date field to R date format
data_video$Date = as.character(data_video$Date)
x <- paste(data_video$Date, data_video$Time_UTC)
date_video_formatR <- strptime(x, "%d/%m/%Y %HH %MM %OS",tz='UTC')
x[which(is.na(date_video_formatR))]
delay <- 15/3*3600/1852;# 15m 3noeud conversion m/s
new_date_video_formatR = date_video_formatR - seconds(delay)

data_EK80_path <- "C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\FAROFA_LatLonBotTime_123FM.csv"    
data_EK80 <- read.csv(data_EK80_path)
colnames(data_EK80)
# Matlab dates conversion into R date format
date_EK80_formatR = as.POSIXct((data_EK80$Time-719529)*86400,origin="1970-01-01",tz='UTC')
# For each GoPro observation, looking for the closest lon/lat from EK80 file


for (i in 1:length(new_date_video_formatR)) {
  i
  # select data of same date 
  data_EK80_small =  data_EK80[which(as.Date(date_EK80_formatR)==as.Date(new_date_video_formatR[i])),]
  date_EK80_formatR_small = date_EK80_formatR[which(as.Date(date_EK80_formatR)==as.Date(new_date_video_formatR[i]))]
  # search minimal value 
  dt = difftime(new_date_video_formatR[i], date_EK80_formatR_small, units = "secs")
  dt = abs(as.numeric(dt))
  ind = which(dt == min(dt))
  
  ## test 
  #new_date_video_formatR[i]
  #date_EK80_formatR_small[ind]
  
  data_video$Longitude[i] = data_EK80_small$Longitude[ind];
  data_video$Latitude[i] = data_EK80_small$Latitude[ind];
  data_video$Bottom[i] = data_EK80_small$Bottom[ind];
}

data_video_all[which(data_video_all$VideoType=="TowedVideos"),] = data_video;

paht2save <- "C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\01-12-2020\\FAROFA_underwater_video_species_observations.csv"
write.csv(data_video_all,paht2save,row.names = FALSE)

pngMAP_df = get_map(location = c(-32.66, -3.92, -32.34,  -3.76), source = "osm", zoom = 12)   
# Checking observations coordinates
x11(width=16, height=12)
map <- ggmap(pngMAP_df)
map <- map + geom_point(data=data_video_all, aes(Longitude,Latitude),shape=21) 
map <- map + theme_bw()
map <- map + labs(title= 'Video sampling')
map
