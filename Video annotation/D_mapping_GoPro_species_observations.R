  #***********************************************************************
  # mapping_GoPro_species_observations.r
  # ------------------------------------
  # SPECIES DISTRIBUTION FROM VIDEOS DATA COLLECTED DURING FAROFA SURVEY
  # ----------------------
  # INPUTS:
  # - Videos data observations (all_DV.csv). See required format at the end of this program
  # - Acoustic data from EK80 echointegration file (FAROFA2.mat) with Time, Longitude and Latitude fields only
  # OUTPUT file :
  # - Videos data observations with coordinates (all_DV_with_coordinates.csv)
  # OUTPUT figures:
  # - GoPro tracks
  # - Species distribution map
  # - Sediment map
  # - Algae map
  #***********************************************************************
  # Creation : Jeremie HABASQUE (IRD)
  # Date : 03/04/2019
  #***********************************************************************
  rm(list=ls())
  
  # Install package 
  install.packages("ggplot2")
  install.packages("mapdata")
  install.packages("scales")
  install.packages("marmap")
  install.packages("reshape2")
  install.packages("ggmap")
  install.packages("R.matlab")
  install.packages("lubridate")
  install.packages("dplyr")
  
  
  # loading librairies
  library(ggplot2)
  library(mapdata)
  library(scales)
  library(marmap)
  library(reshape2)
  library(ggmap)
  library(R.matlab)
  library(lubridate)
  library(dplyr)
  
  #graph.path = 'C:/Users/jhabasqu/Desktop/video_FAROFA_2/'
  
  #reading species observations from video file
  #data_video <- read.csv("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\Big_table_FAROFA_1.csv")
  data_video <- read.csv("C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\Big_table_FAROFA_3.csv")
  
  #head(data_video)
  #str(data_video)
  #summary(data_video)
   
  # convert date field to R date format
  data_video$Date = as.character(data_video$Date)
  data_video$Time_UTC <- seconds_to_period(data_video$Time_UTC)
  x <- paste(data_video$Date, data_video$Time_UTC)
  
  date_video_formatR <- strptime(x, "%Y-%m-%d %HH %MM %OS",tz='UTC')
  x[which(is.na(date_video_formatR))]
  data_video$VideoNumber[which(is.na(date_video_formatR))]
  tmp<-data_video[which(is.na(date_video_formatR)),]
  
  data_video$VideoNumber[which(is.na(date_video_formatR))]
  
  # reading EK80 echointegration file  
  data_EK80 <- readMat('C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie/FAROFA2.mat')
  
  # deleting wrong positions (occurs sometimes with GPS)
  #ind=which(data_EK80$Longitude < -33 | data_EK80$Longitude > -32.3 | data_EK80$Latitude > -3.9 | data_EK80$Latitude < -4);
  #data_EK80$Longitude = data_EK80$Longitude[-ind]
  #data_EK80$Latitude = data_EK80$Latitude[-ind]
  #plot(data_EK80$Longitude, data_EK80$Latitude)
  #data_EK80$Time = data_EK80$Time[-ind];
  #data_EK80$DateTime = data_EK80$DateTime[-ind]
  
  # Matlab dates conversion into R date format
  date_EK80_formatR = as.POSIXct((data_EK80$Time-719529)*86400,origin="1970-01-01",tz='UTC')
  
  # For each GoPro observation, looking for the closest lon/lat from EK80 file
  data_video$Longitude <- NA
  data_video$Latitude <- NA
  data_video$Bottom <- NA
  for (i in 1:length(date_video_formatR)) {
    i
     ind = min(which(abs(as.numeric(difftime(date_video_formatR[i], date_EK80_formatR, units = "secs")))<seconds(1)))
     data_video$Longitude[i] = data_EK80$Longitude[ind];
     data_video$Latitude[i] = data_EK80$Latitude[ind];
     data_video$Bottom[i] = data_EK80$Bottom[ind];
  }
  
  plot(data_video$Longitude, data_video$Latitude)
  data_video$Time <- seconds_to_period(data_video$Time)
  data_video$VideoNumber[which(data_video$Longitude>-32.38 & data_video$Latitude>-3.8)]
  
  # Write new file with coordinates
  #write.csv(data_video,"C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\Big_table_FAROFA_1_coordinates.csv",row.names = FALSE)
  write.csv2(data_video,"C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\FAROFA_underwater_video_species_observations.csv",row.names = FALSE)
  
  data_video <- read.csv("C:\\Users\\jsalveta\\Desktop\\coriger\\Planilhas ajeitadas\\final\\FAROFA_underwater_video_species_observations.csv", sep=",")
  ## MAPS
  
  #coast_map <- fortify(map("worldHires",xlim=c(-33, -32), ylim=c(-4,-3),fill=TRUE, plot=FALSE))
  ## uncomment the following line if you have a Internet connection
  pngMAP_df = get_map(location = c(-32.66, -3.92, -32.34,  -3.76), source = "osm", zoom = 12)   
  #load('C:/Users/garyr/Desktop/jeremy/FdN_map.Rdata')

  # Checking observations coordinates
  x11(width=16, height=12)
  map <- ggmap(pngMAP_df)
  map <- map + geom_point(data=data_video, aes(Longitude,Latitude),shape=21) 
  map <- map + theme_bw()
  map <- map + labs(title= 'Video sampling')
  map
  
  # p <- ggplot(data_video, aes(Longitude,Latitude))
  # p <- p + geom_map(data=coast_map, map=coast_map, aes(x=long, y=lat, map_id=region),fill="grey")
  # p <- p + geom_point() + theme_bw() 
  # p <- p + labs(title= 'Towed GoPro during FAROFA 2')
  # p
  
  # Change name
  
  Config2ChangeColumnNames_path<-"C:\\Users\\jsalveta\\Desktop\\processing video\\list_sp_code.csv"
  Config2ChangeColumnNames <- read.csv(Config2ChangeColumnNames_path,header = TRUE,sep=",")
  newnames <- as.vector(unlist(Config2ChangeColumnNames[,1]))
  oldnames <- as.vector(unlist(Config2ChangeColumnNames[,2]))
  # pull out the names that are actually in x
  old_nms <- oldnames[oldnames %in% names(data_video)]
  new_nms <- newnames[oldnames %in% names(data_video)]
  # rename
  data_video <-  data_video %>% rename_at(vars(old_nms), ~new_nms)
  
  # conversion to a molten data frame 
  data_melt <- melt(data_video,id.vars=c(1:15,47,49,59,68:69,71,73:76), measure.vars = c(16:46,48,50:58,60:67,70,72))#ncol(data_video)-1)
  names<-colnames(data_video)
  colnames(data_melt) <- c(names[1:15],names[47],names[49],names[59],names[68:69],names[71],names[73:76],"Species","Number")
  data_melt <- subset(data_melt, !is.na(Number))
  data_melt$Time <- seconds_to_period(data_melt$Time)
  
  write.csv(data_melt,"C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\\\Melt_table_FAROFA_1_coordinates.csv",row.names = FALSE)
  
  
  data_melt<-read.csv("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas\\bigTables_coordinates\\Melt_table_coordinates.csv")
  
  # x11(width=16, height=12)
  # p <- ggplot(data_melt, aes(Longitude,Latitude))
  # p <- p + geom_map(data=coast_map, map=coast_map, aes(x=long, y=lat, map_id=region),fill="grey")
  # p <- p + geom_point(aes(size=Number)) + scale_size(guide = "none" )
  # p <- p + geom_point(aes(fill=Number,size=Number),shape=21)
  # p <- p +  theme_bw() + facet_wrap(~Species,nrow=3)
  # p <- p + labs(title= 'Species observed by towed GoPro during FAROFA 2')
  # p
  
  x11(width=16, height=12)
  map <- ggmap(pngMAP_df)
  map <- map + geom_point(data=data_melt, aes(Longitude,Latitude, fill=Number, size=Number),shape=21) + scale_size(guide = "none" )
  map <- map + geom_point(data=data_melt, aes(Longitude,Latitude, fill=Species, colour = Species, size=Number,shape = Species)) + scale_size(guide = "none" )
    map <- map +  theme_bw() +  facet_wrap(~Species,nrow=5)
  #map <- map + labs(title= 'Species observed by vertical video profile during FAROFA 1')
  map
  
  x11(width=16, height=12)
  map <- ggmap(pngMAP_df)
  map <- map + geom_point(data=data_melt, aes(Longitude, Latitude, colour=Species,size=Number)) + scale_size(guide = "none" )
  map <- map +  theme_bw() 
  map <- map + labs(title= 'Species observed by vertical video profile during FAROFA 1')
  map

  ############## Sediments map
  
  data_sediments <- subset(data_video, Sediment !="")
  
  x11(width=20, height=15)
  #p <- ggplot(data_sediments, aes(Longitude,Latitude))
  #p <- p + geom_map(data=coast_map, map=coast_map, aes(x=long, y=lat, map_id=region),fill="grey")
  p <- ggmap(pngMAP_df)
  p <- p + geom_point(data=data_sediments, aes(Longitude, Latitude, colour=Sediment)) + scale_size(guide = "none" )
  p <- p + theme_bw() 
  p <- p + labs(title= 'Sediments observed by vertical profile during FAROFA 1')
  p
  
  ############# Algae map 
  
  data_seaweed =subset(data_video, seaweed !="")
  x11(width=16, height=12)
  #p <- ggplot(data_seaweed, aes(Longitude,Latitude))
  #p <- p + geom_map(data=coast_map, map=coast_map, aes(x=long, y=lat, map_id=region),fill="grey")
  p <- ggmap(pngMAP_df)
  p <- p + geom_point(data=data_seaweed, aes(Longitude, Latitude, colour=seaweed)) + scale_size(guide = "none" )
  p <- p + theme_bw() 
  p <- p + labs(title= 'Seaweed observed by towed GoPro during FAROFA 2')
  p
  
  ## REQUIRED FORMAT TO USE THIS PROGRAM 
  
  # [1] "Survey"                 "VideoType"              "number"                 "camera"                
  # [5] "Date"                   "Time"                   "HOUR.UTC"               "HOUR.GOPRO"            
  # [9] "Depth.logbook."         "situation"              "Canthidermis.sufflamen" "Aluterus.scriptus"     
  # [13] "Melichthys.niger"       "unknown"                "Caranx.lugubris"        "Sphyraena.barracuda"   
  # [17] "gelatinous"             "sediment"               "Lactophrys.trigonus"    "shark"                 
  # [21] "Lutjanus.jocu"          "Seriola.dumerili"       "trachurus.sp"           "Caranx.crysos"         
  # [25] "Kyphosus.sectatrix"     "Caranx.latus"           "Acanthurus.coeruleus"   "seaweed"               
  # [29] "Abudefduf.saxatilis"    "Chromis.multilineata"   "Spotlight.parrotfish"   "Aetobatus.narinari"    
  # [33] "Elagatis.bipinnulata"   "Balistes.vetula"        "Comments"  
  