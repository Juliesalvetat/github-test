## Création du tableau pour faire des cartes camemberts ##

# Install package
install.packages("tidyr")
install.packages("dplyr")
install.packages("leaflet")
install.packages("leaflet.minicharts")
install.packages("shiny")
install.packages("htmlwidgets")
install.packages("ggmap")
install.packages("scatterpie")

#Répertoire de travail
setwd("C:\\Users\\jsalveta\\Desktop\\processing video\\Cartography\\Cartography\\Data_Sample")


## CREATION du tableau des max N
Big_table_FAROFA <- read.csv("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas\\Big_table_FAROFA_1.csv", sep=",")

#extraction du tableau contenant uniquement les colonnes avec des espèces dénombrables
especes <- subset(Big_table_FAROFA, select = c(16:46,48,50))
#extraction du vecteur des noms de toutes les espèces 
col<-colnames(especes)

#(Levels renvoie chaque valeur que peut prendre un élement d'une colonne)
#ici on crée un vecteur avec tous les noms possibles des vidéos (VP1, VP2, VP3...)
vp <- levels(Big_table_FAROFA[,4])

colMax <- function(data) sapply(data, max, na.rm = TRUE)

#A CHANGER SI NOMBRE ESPECES CHANGE
table_maxn <- data.frame(matrix(NA, nrow = length(vp), ncol = length(col)))
# On donne leurs noms aux lignes et colonnes
colnames(table_maxn)<-col
rownames(table_maxn)<-vp

# On fait tourner une boucle pour chaque vidéo
for (i in 1:length(vp)) {
  #On utilise le vecteur vp donnant le nom des vidéos pour créer un tableau intermédiaire (b) ne contenant les données que pour une vidéo 
  a = vp[i]
  b <- subset(Big_table_FAROFA, Big_table_FAROFA[,4]==a)
  # On sélectionne uniquementles colonnes des espèces
  c <- subset(b, select = c(16:46,48,50))
  #On applique la fonction maximum
  max <- colMax(c)
  # On transforme le vecteur contenant un nom+une valeur de max dans chaque élement en data frame 
  max <- data.frame(as.list(max))
  # On extrait la ligne contenant les maximums
  max <-as.vector(unlist(max[1,]))
  # On ajoute cette ligne de maximum à la colonne correspondant à la vidéo dans le tableau final 
  table_maxn[i,] <- max # or whatever value you want here
}


### Création du tableau final

# Création d'un tableau vide (1 ligne = 1 vidéo / nb col = Lat Long ID Factor + nb especes)
d <- data.frame(matrix(NA, nrow = length(vp), ncol = 4+length(col)))
colnames(d)<-c("Latitude", "Longitude", "ID", "Factor", col)

# On remplit la colonne ID (correspond au numéro des vidéos) et Factor (correspond à la campagne)
d$ID <- seq(1,length(vp))
d$Factor <- rep("F1")

##Remplissage des colonnes coordonnées

#Importation du tableau des coordonnées
Melt_table_FAROFA_1_coordinates <- read.csv("C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas/Melt_table_FAROFA_1_coordinates.csv")

# Obtention des moyennes des latitudes
vp <- levels(Big_table_FAROFA[,4])

vector_moy_lat = rep(0,length = length(vp))
for (j in 1:length(vp)) {
  #On utilise le vecteur vp donnant le nom des vidéos pour créer un tableau intermédiaire (b) ne contenant les données que pour une vidéo 
  a = vp[j]
  b <- subset(Melt_table_FAROFA_1_coordinates, Melt_table_FAROFA_1_coordinates[,4]==a)
  c = 0
  # Somme des toutes les longitudes
  for (k in 1:length(b)){
    c = c + b[k,18]
  }
  moy_lat = c/length(b)
  vector_moy_lat[j] <- moy_lat
}

#Remplissage colonne latitude    
d$Latitude <- vector_moy_lat

#Idem longitude
vector_moy_long = rep(0,length = length(vp))
for (j in 1:length(vp)) {
  #On utilise le vecteur vp donnant le nom des vidéos pour créer un tableau intermédiaire (b) ne contenant les données que pour une vidéo 
  a = vp[j]
  b <- subset(Melt_table_FAROFA_1_coordinates, Melt_table_FAROFA_1_coordinates[,4]==a)
  c = 0
  # Somme des toutes les longitudes
  for (k in 1:length(b)){
    c = c + b[k,17]
  }
  moy_long = c/length(b)
  vector_moy_long[j] <- moy_long
}

#Remplissage colonne longitude   
d$Longitude <- vector_moy_long

## Remplissage des données pourcentages 

# On créer un tableau vecteur pour y ranger les sommes des maximums par vidéo

#On crée une matrice vide d'une ligne qui contiendra le tableau final 
sum_max <- data.frame(matrix(NA, nrow = 1, ncol = length(vp)))
# On nomme les colonnes et la ligne 
colnames(sum_max)<-vp
rownames(sum_max)<-"Somme_max"

#On fait tourner la boucle pour chaque vidéo 
for (i in 1:length(vp)) {
  #On extrait 1 colonne du tableau des max n ce qui correspond aux données pour une vidéo 
  inter<-as.numeric(table_maxn[i,])
  #is.infinite(inter) > renvoie un vecteur de FALSE et des TRUES s'il y a un maximum
  #inter(is.infinite(inter)) renvoie un vecteur des valeurs correspondantes aux emplacements des true seulement
  # on applique la somme à la liste des max
  sum_max[1,i]=sum(inter[is.finite(inter)])
}

#On créer un tableau pourcent qui contiendra tous les pourcentages d'apparition des maxN 

#Création tableau vide
table_pourcent <- data.frame(matrix(NA, nrow = length(vp), ncol =length(col)))
colnames(table_pourcent)<-especes
rownames(table_pourcent)<-vp

#Remplissage des maximums 
for (i in 1:length(vp)){
  for (j in 1:length(col)){
    table_pourcent[i,j] <- (table_maxn[i,j]/sum_max[,i])
  }
} 

#FINALEMENT : On rentre les pourcentages dans le tableau final 
for (i in 1:length(vp)){
  for (j in 1:length(col)){
    d[i,4+j] <- table_pourcent[i,j]
  }
} 

write.csv(d, file = "table_camembert.csv")

                
