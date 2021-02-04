library(RGeostats)
library(maps)
library(mapdata)

setwd("C:\\Users\\jsalveta\\Desktop\\geostat")
my.dataframe <- read.table(file = 'C:\\Users\\jsalveta\\Desktop\\geostat\\sc_2008_2010.txt')
plot(my.dataframe$Lon,my.dataframe$Lat)
map(add=T)

names(my.dataframe)
my.db <- db.create(my.dataframe)
my.db
my.db <- db.add(my.db,density=Commercial.weight*1000/Area_swept)
my.db
my.db <- db.locate(my.db,c(4,3),loctype="x")
my.db
plot(my.db)
map('worldHires',add=T) #high resolution map

my.db.2010 <- db.sel(my.db,Year==2010) #selection
plot(my.db.2010)
map('worldHires',add=T) #high resolution map


#polygone de travai
poly.sc <- read.table(file = '/amuhome/t18007955/Bureau/a suprimer julie/file03_Polygon_gulf.txt',header=F,sep="\t")
my.poly <-polygon.create(x=-poly.sc[,2],y=poly.sc[,1])
  
?vario.calc
args(vario.calc) 

vg.2010 <- vario.calc(my.db.2010)
plot(vg.2010,npairdw=T,inches=0.05)
vg.2010 <- vario.calc(my.db.2010,lag=0.1,nlag=40,dirvect = c(0,45,90,135))
x11()
plot(vg.2010,npairdw=T,inches=0.05)

vg.2010 <- vario.calc(my.db.2010,lag=0.1,nlag=60)
plot(vg.2010,npairdw=T,inches=0.05)

vg.2010 <- vario.calc(my.db.2010,lag=0.1,nlag=20)
vg.2010.mod.1 <- model.auto(vario=vg.2010,struct = c(1,3,2,12)) #ordre iportant cherche meilleur combinaison
vg.2010.mod.2 <- model.auto(vario=vg.2010,struct = c(1,3,2,12),wmode=0) #wmode=0 fort ponderation

grille <- db.create(flag.grid = T,x0=c(-66,46)+runif(2,-0.0001,0.0001),
                    dx=c(0.025,0.025),nx=c(12,7)*20)# microbruit pour eviter la coincidence entre les pts de grille et les valeurs obs
plot(my.db.2010)
map('worldHires',add=T) #high resolution map
my.polygone <- polygon.digit(col=2)
plot(grille,add=T,fg.in="red")


grille<-db.polygon(grille,my.poly)
plot(grille,add=T,fg.in=2)

#vosinage
args(neigh.create)
voisin <- neigh.create(nmaxi = 20)#valeur par defaut conviennent

kri.2010.1 <- kriging(my.db.2010,grille,vg.2010.mod.1,voisin) # effet de pepite carte plus lisse model precautionneux manque observation fine echelle par rapport a la donnee
kri.2010.2 <- kriging(my.db.2010,grille,vg.2010.mod.2,voisin)

x11()
plot(kri.2010.1,col=rainbow(100,start=0.2,end=1))
map('worldHires',add=T) #high resolution map

x11()
plot(kri.2010.1,name.image=6)#,col=rainbow(100,start=0.2,end=1))
legend.image(zlim = range(kri.2010.1[,6][kri.2010.1[,6]>0]))
plot(my.db.2010,add=T)
map('worldHires',add=T)

x11()
plot(kri.2010.2,col=rainbow(100,start=0.2,end=1))
map('worldHires',add=T) #high resolution map

x11()
plot(kri.2010.2,name.image=6)#,col=rainbow(100,start=0.2,end=1))
legend.image(zlim = range(kri.2010.2[,6][kri.2010.2[,6]>0]))
plot(my.db.2010,add=T)
map('worldHires',add=T)
legend.image(zlim = range(kri.2010.1[,6][kri.2010.1[,6]>0]))

range(kri.2010.1[,6][kri.2010.1[,6]>0]) #intervalle de fluctuation des valeur d'ecart type
range(kri.2010.2[,6][kri.2010.2[,6]>0])
