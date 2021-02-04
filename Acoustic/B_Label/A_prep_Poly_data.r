library(corrplot)


my.dataframe<-read.csv("F123_PolyTable_complete.csv",header = T)
str(my.dataframe)
colnames(my.dataframe)
summary(my.dataframe)
# sardine = my.dataframe[my.dataframe$Poly_LabelName=="small.pelagics.school",]
# summary(sardine)

pngMAP_df = get_map(location = c(-32.66, -3.92, -32.34,  -3.76), source = "osm", zoom = 12)   

x11()
p <- ggplot(my.dataframe,aes(x=Poly_LabelName, fill=FAROFA),) +
  geom_bar(position=position_dodge())+
  #facet_wrap(~FAROFA)+# divide into 2 panels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# 
x11(width=16, height=12)
map <- ggmap(pngMAP_df)
map <- map + geom_point(data=my.dataframe, aes(Poly_Longitude,Poly_Latitude, fill=FAROFA, colour = FAROFA,size=Poly_sum_logSa_70)) + scale_size(guide = "none" )
map <- map +  theme_bw() #+  facet_wrap(~FAROFA)
map <- map + labs(title= 'Small pelagics school')
map

library(ggmap)
pngMAP_df = get_map(location = c(-32.5, -3.92, -32.34,  -3.76), source = "osm", zoom = 12)   
x11(width=16, height=12)
map <- ggmap(pngMAP_df)
map <- map + geom_point(data=sardine, aes(Poly_Longitude,Poly_Latitude),shape=21) 
map <- map + theme_bw()
map <- map + labs(title= 'Small pelagics school')
map


# Supprimer les Label non pertinents
my.dataframe = my.dataframe[!(my.dataframe$Poly_LabelName=="surface.noise" |
                                my.dataframe$Poly_LabelName=="Bottom.corr" |
                                my.dataframe$Poly_LabelName=="correction" |
                                my.dataframe$Poly_LabelName=="Instrument" |
                                my.dataframe$Poly_LabelName=="unknown"|
                                my.dataframe$Poly_LabelName=="bottom.high200"),]
str(my.dataframe)
my.dataframe$FAROFA = as.factor(my.dataframe$FAROFA)
my.dataframe$MPA= as.factor(my.dataframe$MPA)
my.dataframe$Wind_exposure = as.factor(my.dataframe$Wind_exposure)

my.dataframe$Poly_LabelName <- factor(my.dataframe$Poly_LabelName)

#Passage au log de sum Sa
my.dataframe$Poly_sum_logSa_70 = log10(my.dataframe$Poly_sum_Sa_70+1)
my.dataframe$Poly_sum_logSa_200 = log10(my.dataframe$Poly_sum_Sa_200+1)

my.dataframe$logSaFish_70 = log10(my.dataframe$SaFish_70+1)
my.dataframe$logSaFish_200 = log10(my.dataframe$SaFish_200+1)

my.dataframe$logSaNoFish_70 = log10(my.dataframe$SaNoFish_70+1)
my.dataframe$logSaNoFish_200 = log10(my.dataframe$SaNoFish_200+1)

# Ajouter Sa Fish other 

my.dataframe$SaFish_Other_70 = abs(my.dataframe$SaFish_70 - my.dataframe$Poly_sum_Sa_70)
my.dataframe$SaFish_Other_200 = abs(my.dataframe$SaFish_200 - my.dataframe$Poly_sum_Sa_200)

my.dataframe$logSaFish_Other_70 = log10(my.dataframe$SaFish_Other_70+1)
my.dataframe$logSaFish_Other_200 = log10(my.dataframe$SaFish_Other_200+1)


#correction distance au shelf_break

my.dataframe$distanceToShelfbreak_factor = my.dataframe$distanceToShelfbreak
my.dataframe$distanceToShelfbreak_factor[my.dataframe$Poly_Bottom <=20]="eu_upper"
my.dataframe$distanceToShelfbreak_factor[my.dataframe$Poly_Bottom >20 & my.dataframe$Poly_Bottom <=40]="eu_lower"
my.dataframe$distanceToShelfbreak_factor[my.dataframe$Poly_Bottom >40 & my.dataframe$Poly_Bottom <=60]="meso_upper"
my.dataframe$distanceToShelfbreak_factor[my.dataframe$Poly_Bottom >60 & my.dataframe$Poly_Bottom <=80]="meso_middle"
my.dataframe$distanceToShelfbreak_factor[my.dataframe$Poly_Bottom >80]="meso_lower"

my.dataframe$distanceToShelfbreak_factor = as.factor(my.dataframe$distanceToShelfbreak_factor)
str(my.dataframe)

my.dataframe$distanceToShelfbreak[my.dataframe$Poly_Bottom <50 ]=-abs(my.dataframe$distanceToShelfbreak[my.dataframe$Poly_Bottom <50 ])

#rename LabelNames
rename.data.frame = read.csv("rename-label.csv",header = T)
bug = cbind(levels(my.dataframe$Poly_LabelName), as.vector(unlist(rename.data.frame[,2])))
#my.dataframe = na.omit(my.dataframe)
levels(my.dataframe$Poly_LabelName)<-as.vector(unlist(rename.data.frame[,2]))
#my.dataframe = na.omit(my.dataframe)
my.dataframe$Poly_LabelName <- factor(my.dataframe$Poly_LabelName)

## enleve la drina 
my.dataframe = my.dataframe[my.dataframe$Poly_Longitude >-32.55,]; # enleve la Drina 
my.dataframe$Poly_LabelName <- factor(my.dataframe$Poly_LabelName)

## enleve premiere column
my.dataframe = my.dataframe[,-c(1)]

write.csv(my.dataframe,'./modified data/my_Poly_data.csv',row.names = F)

######

my.dataframe = read.csv('./modified data/my_Poly_data.csv',header = T,stringsAsFactors=T)
my.dataframe$FAROFA = as.factor(my.dataframe$FAROFA)
str(my.dataframe)
colnames(my.dataframe)

my.dataframe_reduced = my.dataframe[,(colnames(my.dataframe)=="Poly_Longitude"|
                                         colnames(my.dataframe)=="Poly_Latitude"|
                                         colnames(my.dataframe)=="Poly_LabelName"|
                                         colnames(my.dataframe)=="Poly_sum_logSa_70"|
                                        colnames(my.dataframe)=="Poly_sum_logSa_200"|
                                         colnames(my.dataframe)=="logSaFish_70"|
                                         colnames(my.dataframe)=="logSaNoFish_70"|
                                        colnames(my.dataframe)=="logSaFish_200"|
                                        colnames(my.dataframe)=="logSaNoFish_200"|
                                        colnames(my.dataframe)=="logSaFish_Other_70"|
                                        colnames(my.dataframe)=="logSaFish_Other_200"|
                                         colnames(my.dataframe)=="Poly_Bottom"|
                                         colnames(my.dataframe)=="Rugosity"|
                                         colnames(my.dataframe)=="Slope"|
                                         colnames(my.dataframe)=="distanceToCoast"|
                                         colnames(my.dataframe)=="distanceToShelfbreak"|                                        
                                         colnames(my.dataframe)=="FAROFA"|
                                         colnames(my.dataframe)=="MPA"|
                                         colnames(my.dataframe)=="Wind_exposure"|
                                         colnames(my.dataframe)=="Sediment"|
                                         colnames(my.dataframe)=="distanceToShelfbreak_factor")]                                                 





# str(my.dataframe_reduced)
# res<-cor(my.dataframe_reduced)
# col<-colnames(my.dataframe_reduced)
# x11()
# corrplot(res, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)

my.dataframe_reduced = my.dataframe_reduced[,(colnames(my.dataframe_reduced)=="Poly_sum_logSa_70"|
                                              colnames(my.dataframe_reduced)=="logSaFish_70"|
                                              colnames(my.dataframe_reduced)=="logSaNoFish_70"|
                                                colnames(my.dataframe_reduced)=="logSaFish_Other_70"|
                                                colnames(my.dataframe_reduced)=="logSaFish_Other_200"|
                                                colnames(my.dataframe_reduced)=="logSaNoFish_200"|
                                              colnames(my.dataframe_reduced)=="Poly_Bottom"|
                                              colnames(my.dataframe_reduced)=="Rugosity"|
                                              colnames(my.dataframe_reduced)=="Slope"|
                                              colnames(my.dataframe_reduced)=="distanceToCoast"|
                                              colnames(my.dataframe_reduced)=="distanceToShelfbreak")] 

my.dataframe_quali = my.dataframe[,(colnames(my.dataframe)=="Poly_LabelName"|
                                      colnames(my.dataframe)=="FAROFA"|
                                      colnames(my.dataframe)=="MPA"|
                                      colnames(my.dataframe)=="Wind_exposure"|
                                      colnames(my.dataframe)=="Sediment"|
                                      colnames(my.dataframe)=="distanceToShelfbreak_factor")]  

str(my.dataframe_reduced)
str(my.dataframe_quali)

res<-cor(my.dataframe_reduced)
corrplot(res, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45)


pairs(my.dataframe_reduced)

# 
# 
# heatmaply(
#   my.dataframe_reduced_pca,
#   xlab = "Features",
#   ylab = "Label", 
#   main = "raw data"
# )

my.dataframe_reduced_pca = cbind(my.dataframe_quali,my.dataframe_reduced)
str(my.dataframe_reduced_pca)
write.csv(my.dataframe_reduced_pca,'./modified data/my_Poly_data_pca.csv',row.names = F)



###################### modify ESU 
my.dataframe_ESU<-read.csv("F123_Averaged_SaSum_Sed_complete.csv",header = T)
str(my.dataframe_ESU)
colnames(my.dataframe_ESU)

# my.dataframe = my.dataframe_ESU[!is.na(my.dataframe_ESU$Sediment),]
# my.dataframe = my.dataframe[which(my.dataframe$Sediment!="Un"),]
# 
# x11(width=16, height=12)
# map <- ggmap(pngMAP_df)
# map <- map + geom_point(data=my.dataframe, aes(Lon_mean,Lat_mean, fill=Sediment, colour = Sediment)) + scale_size(guide = "none" )
# map <- map +  theme_bw() #+  facet_wrap(~FAROFA)
# map <- map + labs(title= 'Sediment')
# map



## enleve la drina 
my.dataframe_ESU = my.dataframe_ESU[my.dataframe_ESU$Lon_mean >-32.55,]; # enleve la Drina 

my.dataframe_ESU$logSaFish_70 = log10(my.dataframe_ESU$SaFish_70+1)
my.dataframe_ESU$logSaFish_200 = log10(my.dataframe_ESU$SaFish_200+1)
my.dataframe_ESU$logSaNoFish_70 = log10(my.dataframe_ESU$SaNoFish_70+1)
my.dataframe_ESU$logSaNoFish_200 = log10(my.dataframe_ESU$SaNoFish_200+1)
my.dataframe_ESU$distanceToShelfbreak_factor = my.dataframe_ESU$distanceToShelfbreak
my.dataframe_ESU$distanceToShelfbreak_factor[my.dataframe_ESU$Poly_Bottom <=20]="eu_upper"
my.dataframe_ESU$distanceToShelfbreak_factor[my.dataframe_ESU$Poly_Bottom >20 & my.dataframe_ESU$Poly_Bottom <=40]="eu_lower"
my.dataframe_ESU$distanceToShelfbreak_factor[my.dataframe_ESU$Poly_Bottom >40 & my.dataframe_ESU$Poly_Bottom <=60]="meso_upper"
my.dataframe_ESU$distanceToShelfbreak_factor[my.dataframe_ESU$Poly_Bottom >60 & my.dataframe_ESU$Poly_Bottom <=80]="meso_middle"
my.dataframe_ESU$distanceToShelfbreak_factor[my.dataframe_ESU$Poly_Bottom >80]="meso_lower"

my.dataframe_ESU$distanceToShelfbreak_factor = as.factor(my.dataframe_ESU$distanceToShelfbreak_factor)
str(my.dataframe_ESU)

my.dataframe_ESU$distanceToShelfbreak[my.dataframe_ESU$Poly_Bottom <50 ]=-abs(my.dataframe_ESU$distanceToShelfbreak[my.dataframe_ESU$Poly_Bottom <50 ])


# my.dataframe_ESU = my.dataframe_ESU[,(colnames(my.dataframe_ESU)=="Depth_mean"|
#                                                 colnames(my.dataframe_ESU)=="logSaNoFish_70"|
#                                                 colnames(my.dataframe_ESU)=="logSaFish_70"|
#                                                 colnames(my.dataframe_ESU)=="logSaNoFish_200"|
#                                                 colnames(my.dataframe_ESU)=="logSaFish_200"|
#                                                 colnames(my.dataframe_ESU)=="Rugosity"|
#                                                 colnames(my.dataframe_ESU)=="Slope"|
#                                                 colnames(my.dataframe_ESU)=="distanceToCoast"|
#                                                 colnames(my.dataframe_ESU)=="distanceToShelfbreak")] 

#pairs(my.dataframe_ESU)
write.csv(my.dataframe_ESU,'./modified data/my_ESU_table.csv',row.names = F)


test3<-head(my.dataframe_ESU[order(my.dataframe_ESU$Depth_std,decreasing=TRUE),])
head(my.dataframe_ESU[order(my.dataframe_ESU$Slope,decreasing=TRUE),])
head(my.dataframe_ESU[order(my.dataframe_ESU$Rugosity,decreasing=TRUE),])

my.dataframe_ESU$Depth_diff = abs(my.dataframe_ESU$Depth_start-my.dataframe_ESU$Depth_end)
head(my.dataframe_ESU[order(my.dataframe_ESU$Depth_diff,decreasing=TRUE),])

test<-my.dataframe_ESU[c(23590:23600),]
max(my.dataframe_ESU$Depth_diff)
