library(RGeostats)
library(maps)
library(mapdata)
library(RColorBrewer) ; my.palette <- brewer.pal(9,"YlOrRd")

source("C:/Users/jsalveta/Documents/1_WORK/Cours/geostat/transformation_des_coordonnees_par_trapeze 2020_02_07.txt")

my.dataframe <- read.csv(file = 'Depth_Charts.csv', sep=';')

#regroupe les 2 ABRACOS et les 3 FAROFA
my.dataframe$col = my.dataframe$Source
levels(my.dataframe$col) <- c("Charts","ABRACOS","ABRACOS","FAROFA","FAROFA","FAROFA","Lines")  
#windows()
#plot(my.dataframe$lon_shape,my.dataframe$lat_shape, col=my.dataframe$col)
#legend('topleft',legend = levels(my.dataframe$col),col=1:4,pch=1)

#my.dataframe = my.dataframe [-which(my.dataframe$Source=="Charts"|my.dataframe$Source=="Lines"),]
my.dataframe = my.dataframe [which(my.dataframe$col=="FAROFA"),] # considere seulement les donnees issues de FAROFA 
my.dataframe = my.dataframe[my.dataframe$lon_shape >-32.55,]; # enleve la Drina 
my.dataframe = my.dataframe[my.dataframe$Depth_Points < 60,]; # enleve le talu 

my.dataframe = my.dataframe[my.dataframe$Depth_Points != 50.7998148,]; # enleve 170 valeurs suspectes car ?gales ? 50.7998148 
#                                                                et ? cot? de valeurs bcp plus faibles

tmp <- (my.dataframe$Source=="Farofa1") + (my.dataframe$Source=="Farofa2")*2 + (my.dataframe$Source=="Farofa3")*3
my.dataframe$col <- tmp

# lecture des polygones
my.polygone_shelf <- readRDS("my.polygone_shelf.rds")
tmp <- my.polygone_shelf@sets[[1]]
poly.Fernando <- polygon.create(tmp$x[c(1:25,48,49)],tmp$y[c(1:25,48,49)])
poly.shelf <- polygon.create(tmp$x[26:47],tmp$y[26:47])
# reduction de Fernando ? un perimetre convexe clockwise ferme
tmp <- chull(poly.Fernando@sets[[1]]$x[c(1:24)],poly.Fernando@sets[[1]]$y[c(1:24)])
sel.ref <- rev(tmp[order(tmp)]) # clockwise
sel.ref <- c(sel.ref,sel.ref[1]) # fermeture
ref.line <- list(x=poly.Fernando@sets[[1]]$x[sel.ref],y=poly.Fernando@sets[[1]]$y[sel.ref])

#my.dataframe = na.omit(my.dataframe)
#windows()
#plot(my.dataframe$lon_shape,my.dataframe$lat_shape, col=my.dataframe$Source)
#legend('topleft',legend = levels(my.dataframe$Source),col=1:3,pch=1)

# regularisation des donn?es
my.db <- db.create(my.dataframe)
#my.db <- db.add(my.db,density=Depth_Points)
my.db <- db.locate(my.db,c(3,4),loctype="x")
my.db <- db.locate(my.db,"Depth_Points",loctype="z")
regul.grid <- db.grid.init(my.db,nodes=500)
regul.grid <- db.stat.grid(my.db,regul.grid)
#vario <- vario.grid(regul.grid,nlag=200)
#model <- model.auto(vario,c(8,3),wmode=0)
# ecriture des donn?es regularisees sous forme de point et non de grille
sel <- !is.na(regul.grid[,4])
regul.pt <- db.create(x1=regul.grid[,2][sel],x2=regul.grid[,3][sel],z1=regul.grid[,4][sel])
# d?finition de la grille de krigeage
grid <- db.grid.init(regul.pt,nodes=150)
#grid <- db.selhull(my.db,grid,combine = "AND")
grid <- db.polygon(grid,poly.shelf,combine = "and")
grid <- db.polygon(grid,poly.Fernando,flag.out=T,combine = "and")
# D?finition du voisinage
#kri <- kriging(regul.pt,grid,model,neigh)
#plot(kri)
#plot(regul.pt,add=T,inches=0.5)

# projection des points 
# suppression des points ? l'interieur de l'enveloppe convexe
regul.pt <- db.polygon(regul.pt,polygon.create(ref.line$x,ref.line$y),flag.out = T)
regul.pt <- db.reduce(regul.pt)
regul.conform <- proj.conformal(regul.pt[,2],regul.pt[,3],ref.line,orthogonal = F,closed=T)
db.regul.conform <- db.create(regul.conform)
db.regul.conform <- db.locate(db.regul.conform,2:3,"x")
db.regul.conform <- db.add(db.regul.conform,z1=regul.pt[,4])
plot(db.regul.conform)

# projection de la grille
# extraction des points actifs de la grille
df <-db.extract(grid,2:3,flag.compress = T) 
# projection ATTENTION l'objet resultat s'appelle grid.conform mais n'est plus une grille ? cause de la projection
grid.conform <- proj.conformal(df$x1,df$x2,ref.line,orthogonal = F,closed=T)
db.grid.conform <- db.create(grid.conform)
db.grid.conform <- db.locate(db.grid.conform,2:3,"x")

# vario conform
vario.conform <- vario.calc(db.regul.conform,lag=c(0.0001,0.0025),nlag=c(300,15),dirvect=c(0,90))
#breaks=list(breaks.perp,breaks.along),dirvect = c(0,90))
plot(vario.conform,npairdw=T)

# model conform
model.conform <- model.auto(vario.conform,struct=c(1,3,3,12,12),wmode=3)

#voisinage
#neigh <- neigh.create(ndim=2,nmini=20,nmaxi = 100,radius=0.02)
rayon <- 0.02
neigh.conform <- neigh.create(type=2,radius=rayon,flag.continuous = T,flag.sector = T,nsect = 6,nsmax=40)

# ajout de blocs de donn?es aux deux extr?mit?s pour assurer la continuit? spatiale qu'on a cass?e avec la projection
yrange <- range(db.regul.conform[,3])
#selection des blocs
bloc.bas <- db.reduce(db.sel(db.regul.conform, y < rayon))
bloc.haut <- db.reduce(db.sel(db.regul.conform, y > (yrange[2] - rayon)))
# d?calage avant raccordement : la valeur 0.001 est choisie arbitrairement 
# elle doit permettre de raccrocher de fa?on coh?rente les donn?es .....
bloc.haut[,3] <- bloc.haut[,3]-yrange[2] - 0.001
bloc.haut <- db.extract(bloc.haut,2:4)
bloc.bas[,3] <- bloc.bas[,3]+yrange[2] + 0.001
bloc.bas <- db.extract(bloc.bas,2:4)

db.regul.conform.expanded <- db.append(db.regul.conform,bloc.haut)
db.regul.conform.expanded <- db.append(db.regul.conform.expanded,bloc.bas)


# attention avec 150 noeuds dans la d?finition de la grille ca prend 5 minutes 
kri.conform <- kriging(db.regul.conform.expanded,db.grid.conform,model.conform,neigh.conform)

kri <- grid
sel <- grid[,5]
tmp <- rep(NA,grid$nech)
tmp[sel] <- kri.conform[,4]
kri <- db.add(kri,Kriging = tmp)
plot(kri,col=my.palette)

varkri = kri@items[["Kriging"]]
lon = kri@items[["x1"]]
lat = kri@items[["x2"]]

matkri= matrix(varkri,150,150)  
matlon= matrix(lon,150,150)  
matlat= matrix(lat,150,150)  

table_kri = cbind(lon,lat,varkri)
write.csv(matkri,"matkribathy.csv")
write.csv(matlon,"matlonbathy.csv")
write.csv(matlat,"matlatbathy.csv")

# reste ? faire dans le cadre de cette fa?on de faire (il existe des alternatives eg SPDE, coK des farofa_1,2,3) :
# - choisir un voisinage anisotrope
# - augmenter la d?finition de la grille de krigeage (attention au temps de calcul .... 
# ? ne faire qu'une fois la chaine de traitement bien fix?e)
# jouer avec la representation graphique de la sortie

