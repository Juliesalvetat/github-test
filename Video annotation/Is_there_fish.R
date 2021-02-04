library(janitor)
library(lubridate)

colMax <- function(data) sapply(data, max, na.rm = TRUE)

#path2video <- "C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\VerticalProfiles\\2018-04-19\\VP_06\\Processed_GOPR0482.csv"
#path2video <- "C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\VerticalProfiles\\2018-04-23\\VP_24\\Processed_GOPR0553.csv"
#path2video <- "C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\VerticalProfiles\\2018-04-23\\VP_24\\Processed_GOPR0554.csv"

#path2video <- "C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-19\\DV_13\\Processed_GOPR7448.csv"
#path2video <- "C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-19\\DV_13\\Processed_GP017448.cv"



#path2video <-"C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-20\\DV_14\\Processed_GOPR9575.csv"
#path2video <-"C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-20\\DV_14\\Processed_GP019575.csv"
#path2video <-"C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-20\\DV_14\\Processed_GP029575.csv"
#path2video <-"C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\FAROFA_2\\CSV\\TowedVideos\\2018-04-20\\DV_14\\Processed_GP039575.csv"

video_csv <- read.csv(path2video,header = TRUE, sep = ",", quote = "\"", 
                               dec = ".", fill = FALSE, comment.char="")
list_sp_code<-read.csv("C:\\Users\\stagiaire\\Desktop\\Planilhas ajeitadas\\list_sp_code.csv",header = TRUE, sep = ",", quote = "\"", 
         fill = FALSE, comment.char="")

# convert Time in minute to find fish on video more easily
video_csv$Time<-lubridate::seconds_to_period(video_csv$Time)

test = remove_empty(video_csv, which = c("cols")) # remove empty columns
test2 = remove_empty(test[c(1:6,9:length(test))], which = c("rows")) #select fish columns
df_fish = test2[!rowSums(is.na(test2)) > ncol(test2)-7,] # remove rows with more then ncol(test2)-7 Na, remove rows where there is no fish

c<-data.frame(colMax (df_fish[c(7:length(df_fish))]))