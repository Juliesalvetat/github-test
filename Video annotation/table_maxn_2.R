## TABLEAU DES MAX N ##
# Input : Dossier FAROFA_1 avec tous les dossiers Processe_G....csv 
# Output : Tableau avec le max de chaque espèces pour chaque vidéo 


#On charge le répertoire de travail comme étant FAROFA_1
setwd("/Users/leadoucet/Desktop/Rdonnee")

#Création de la fonction qui retire le max de chaque colonne
colMax <- function(data) sapply(data, max, na.rm = TRUE)

#Importe la table réunissant toutes les observations de poissons pour toutes les vidéos
Big_table_FAROFA_2 <- read.csv("C:\\Users\\jsalveta\\Desktop\\Suelen\\Videos\\Video_annotations_FAROFA_2_with_coordinates.csv", sep=",")
#extraction des colonnes espèces 
especes <- subset(Big_table_FAROFA_2, select = c(18:length(Big_table_FAROFA_2)-3))

#VP1 <- subset(Big_table_FAROFA_1, Big_table_FAROFA_1[,4]=="VP_01")
#Levels renvoie chaque valeur que peut prendre un élement d'une colonne, ici on crée un vecteur avec VP1, vp2, VP3...
vp <- levels(Big_table_FAROFA_2[,4])
             
for (i in range(1:17)) {
  table_maxn <- data.frame(,row.names = )
  a = vp[i]
  b <- subset(Big_table_FAROFA_2, Big_table_FAROFA_2[,4]==a)
  c <- subset(b, select = c(14:length(b)))
  max <- colMax(c)
  table_maxn <- c(table_maxn, max)
}






#on applique la fonction max à process pour obtenir le max de chaque colonne. ca affiche -Inf quand il n'y a pas de l'espèce
maxprocess <- colMax(process)
print(maxprocess)
# création du dataframe pour ranger les espèces
maxF1 <- data.frame(VP1 = maxprocess, row.names=colnames(process))
print(maxF1)