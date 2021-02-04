library("plyr") 
library(RGeostats)
library(maps)
library(mapdata)
library(wesanderson)
# Definition d'une palette
my.palette = wes_palette("Zissou1", 16, type = "continuous")
new_name_and_path<-"C:\\Users\\jsalveta\\Desktop\\Planilhas ajeitadas - Copie\\Planilhas ajeitadas\\FAROFA_melSA.csv"
my.dataframe<-read.csv(new_name_and_path,header = T)
setwd("c:/Users/jsalveta/Documents/1_WORK/Cours/geostat")
source("C:/Users/jsalveta/Documents/1_WORK/Cours/geostat/transformation_des_coordonnees_par_trapeze 2020_02_07.txt")


my.dataframe = na.omit(my.dataframe)
my.dataframe$Horizontal_logSa70 = log10(my.dataframe$Horizontal_sa_surface_fish_70+1)
my.dataframe$Horizontal_logSa200 = log10(my.dataframe$Horizontal_sa_surface_fish_70+1)
my.dataframe = my.dataframe [my.dataframe$Longitude>-33,] 
my.dataframe = my.dataframe [my.dataframe$Latitude < -3.5,]
my.dataframe$time_formatR = as.POSIXct((my.dataframe$Time-719529)*86400,origin="1970-01-01",tz='UTC')
carre_formatR <- strptime("2019-04-18 14:58:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
carre_formatR_fin <- strptime("2019-04-20 10:00:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
my.dataframe = my.dataframe[!(my.dataframe$time_formatR>carre_formatR & my.dataframe$time_formatR<carre_formatR_fin),]

my.dataframe[which(my.dataframe$Horizontal_sa_surface_fish_70>400),]

my.dataframe[which(my.dataframe$Horizontal_sa_surface_fish_70>400000),]

max(my.dataframe$Horizontal_sa_surface_fish_70)
# creation DB
db <- db.create(my.dataframe)
db <- db.locate(db,c(2,3),loctype="x")
db <- db.locate(db,c("Horizontal_sa_surface_fish_70"),"z")

# Simplification des noms des variables
db <- db.rename(db,c("Horizontal_sa_surface_fish_70","Horizontal_logSa70","Bottom"),c("Sa70","LogSa70","Depth"))

# Regularisation
# Une fois acquis le fait que les differentes FAROFA peuvent etre melangees (resultats not shown), on fait des moyennes de Sa par pixel
# en melangeant les 3 FAROFA
pixel.size <- 0.0005
regul.grid <- db.grid.init(db,dcell=pixel.size)#nodes=150)
regul.grid <- db.stat.grid(db,regul.grid,names=c("Sa70","LogSa70","Depth","FAROFA"),radix="",modify.target = F)
regul.grid <- db.locate(regul.grid,"Sa70","z")
regul.grid <- db.rename(regul.grid,c("x1","x2"),c("Longitude","Latitude"))
x11()
plot(regul.grid,col=my.palette,asp=1)

# # Ecriture des donnees regularisees sous forme de point et non de grille.
# # Cela ppermet de faire des calculs de variogrammes avec tolerance sur la direction.
# # vario.grid() ne calcule que le long des axes principaux de la grille et avec tous les points de la grille
# # y compris toutes les valeurs NA !!!
regul.pt <- db.grid2point.sampling(regul.grid,npack=1,names=c("Sa70","LogSa70","Depth","FAROFA"))
regul.pt <- db.reduce(db.sel(regul.pt,!is.na(LogSa70)))
regul.pt <- db.locate(regul.pt,c("LogSa70","Depth","FAROFA"),NA)
x11()
plot(regul.pt,col=my.palette,zmin=0.0001,zmax=300,asp=1)
#my.polygone <- polygon.digit(col=2)
#saveRDS(my.polygone,"my.polygone_baliste.rds")
my.polygone=readRDS("my.polygone_baliste.rds")

x11()
plot(db.polygon(regul.grid,my.polygone),col=my.palette,asp=1)





########### estimation de biomasse presence absence

db.regul.baliste <- db.reduce(db.polygon(regul.grid,my.polygone))
db.baliste <- db.reduce(db.polygon(db,my.polygone))

plot(db.regul.baliste,col=my.palette,asp=1)

plot(db.sel(db.regul.baliste,Sa70==0),asp=1,col=1)
plot(db.sel(db.regul.baliste,Sa70>0),asp=1,add=T)

db.regul.baliste <- db.add(db.regul.baliste,PA=Sa70>0)

db.regul.baliste <- db.locate(db.regul.baliste,5,"z")
vario.baliste <- vario.calc(db.sel(db.regul.baliste,Sa70>0), lag=0.001,nlag=50)
plot(vario.baliste)

db.regul.baliste <- db.locate(db.regul.baliste,9,"z")
vario.PA <- vario.calc(db.regul.baliste, lag=0.001,nlag=50,dirvect=c(0,90))
plot(vario.PA)

vario.mod <- model.auto(vario.PA,struct=c(1,2))

grille <- db.grid.init(db.regul.baliste,nodes=150)
neigh <- neigh.create(nmini=20,nmaxi=50,radius = 0.02)

krigeage.PA <- kriging(db.regul.baliste,db.polygon(grille,my.polygone),vario.mod,neigh)
plot(krigeage.PA,asp=1)

krigeage.PA <- db.add(krigeage.PA,estim.PA=Kriging.PA.estim >0.5)
plot(krigeage.PA,asp=1)
plot(my.polygone,add=T)

surface.presence <- sum(krigeage.PA[,7]==1,na.rm=T)*(grille@dx[1]*grille@dx[2])*60*60*1852^2

tmp <- db.sel(db.regul.baliste,Sa70>0)
tmp2 <- db.extract(tmp,4)
m <- mean(tmp2) #valeur moyenne des Sa positives
abondance <- m*surface.presence/(1852^2*4*pi*10^(-39.3/10)) #en nbr d'individus
poids.moyen <- 100 #gramme
biomasse <- abondance * poids.moyen /10^6


#####################################################3

#######################
#######################
# Passage aux indicatrices
#######################
#######################

# Ajout des indicatrices ? la base de donn?es regularisees
# Ind1 = presence/absence
# Ind2-4 = quartile des donn?es positives
db.regul.baliste <- db.reduce(db.polygon(regul.grid,my.polygone))
db.regul.baliste <- db.reduce(db.sel(db.regul.baliste,!is.na(Sa70)))
db.regul.baliste <- db.delete(db.regul.baliste,10)
db.regul.baliste <- db.locate(db.regul.baliste,4,"z")
plot(db.sel(db.regul.baliste,Sa70==0),asp=1,col=1)
plot(db.sel(db.regul.baliste,Sa70>0),asp=1,add=T)

# dbin <- db.regul.baliste
# zcut <- quantile(dbin$Sa70[dbin$Sa70 > 0],seq(0.25,1,0.25),na.rm=T)
# my.limits <- limits.create(mini=c(0,zcut[-4]),maxi=zcut,flag.zcut.int=flag.zcut.int,incmini=T,incmaxi=c(F,F,F,T))
# print(my.limits)
# dbout <- db.indicator(dbin, limits=my.limits,radix="Ind")
# db.regul.baliste.ind <-  dbout

dbin <- db.regul.baliste
zcut <- quantile(dbin$Sa70[dbin$Sa70 > 0],seq(0.25,1,0.25),na.rm=T)
my.limits <- limits.create(mini=c(0,0.01,zcut[-4]),maxi=c(0.01,zcut),incmini=c(T,rep(F,4)),incmaxi=rep(T,5),flag.zcut.zero=T)
print(my.limits)
dbout <- db.indicator(dbin, limits=my.limits,radix="Ind")
db.regul.baliste.ind <-  dbout

plot(db.sel(db.regul.baliste.ind,Ind.Sa70.1==1),inches=0.1,pch=20,col=my.palette[1],asp=1)
plot(db.sel(db.regul.baliste.ind,Ind.Sa70.2==1),inches=0.1,pch=20,col=my.palette[5],add=T)
plot(db.sel(db.regul.baliste.ind,Ind.Sa70.3==1),inches=0.1,pch=20,col=my.palette[11],add=T)
plot(db.sel(db.regul.baliste.ind,Ind.Sa70.4==1),inches=0.1,pch=20,col=my.palette[16],add=T)

### MAF
# on travaille avec 3 indicatrices ; la 4eme se d?duisant des 3 premi?res.
db.regul.baliste.ind <- db.locate(db.regul.baliste.ind,names=10:14, loctype = "z", flag.locnew = TRUE)
maf <- maf.calc(db.regul.baliste.ind,h0=pixel.size,dh=pixel.size/2)
db.regul.baliste.ind.maf <- pca.z2f(db.regul.baliste.ind,pca=maf,verbose = F)

# variogramme des MAFs
vario.maf <- vario.calc(db.regul.baliste.ind.maf,lag=pixel.size,nlag=15)#,dirvect = c(0,90)) 
x11()
plot(vario.maf,npairdw=T,col=my.palette[c(1,16)])

# extraction des variogrammes simples en vue de faire des Krigeages Simples
for(i in 1:5) assign(paste("vario.maf.",i,sep=""),vario.reduce(vario.maf,varcols=i,flag.vario=T))
for(i in 1:5) plot(get(paste("vario.maf.",i,sep="")),add=i>1,ylim=c(0,1.4),npairdw=T,col=my.palette[c(1,16)])


# modelisation des variogrammes
for(i in 1:5) assign(paste("model.maf.",i,sep=""),model.auto(get(paste("vario.maf.",i,sep="")),struct=c(1,3,3),wmode=2,npairdw=T))

# r?alisation des krigeages monovariables des 3 MAFs*
# adaptation de la taille du voisinage ? la port?e des mod?les
# attention si le voisinage est trop restreint il apparait des zones non estim?es ...
rayon <- c(0.02,0.01,0.01,0.01,0.01)

grille <- db.grid.init(db.regul.baliste,nodes=200)

# temps de clacul des 3 krigeages monovariable environ 10-20 mn.
# en cokrigeage multivari? la nuit ne suffit pas ...
for(i in 1:5){
  cat("kriging MAF ",i,"\n")
  neigh.conform <- neigh.create(type=2,radius=rayon[i],flag.continuous = T,
                                flag.sector = T,nsect = 6,nsmax=40,
                                nmini=10,nmaxi=240,
                                flag.aniso=F)
  db.regul.baliste.ind.maf <- db.locate(db.regul.baliste.ind.maf,paste("Factor.",i,sep=""),"z", flag.locnew = TRUE)
  res <- kriging(db.regul.baliste.ind.maf,db.polygon(grille,my.polygone),model=get(paste("model.maf.",i,sep="")),neigh=neigh.conform)
  assign(paste("kri.maf.",i,sep=""),res)
}

plot(db.regul.baliste.ind.maf)

kri <- kri.maf.1
kri <- db.add(kri,Kriging.Factor.2.estim=kri.maf.2[,5],
              Kriging.Factor.3.estim=kri.maf.3[,5],
              Kriging.Factor.4.estim=kri.maf.4[,5],
              Kriging.Factor.5.estim=kri.maf.5[,5])

plot(kri.maf.2,col=my.palette,pos.legend=3)
plot(poly.Fernando,add=T)

### Reconstruction des densit?s interpol?es ? partir des MAFs interpol?s
kri <- db.locate(kri,c(5,7:10),"z")
res <- pca.f2z(kri,pca=maf,verbose = F)

# seuillage des krigeages d'indicatrices ? [0,1]
for(i in 11:15){
  res[,i][res[,i] > 1] <- 1
  res[,i][res[,i] < 0] <- 0
}

# plot(res,name=9,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(Sa < q25+)",asp=1)
# plot(res,name=10,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(q25+ < Sa < q50+)")
# plot(res,name=11,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(Sa > q75+)",asp=1,las=1) ; plot(poly.Fernando,add=T)

##### a supprimer res <- db.add(res,Variable.4=1-(Variable.1+Variable.2+Variable.3))

res <- db.compare(res,fun="maxi",names=11:15)
res <- db.add(res, classe = 1*(Variable.1==maxi)+2*(Variable.2==maxi)+3*(Variable.3==maxi)+4*(Variable.4==maxi)+5*(Variable.5==maxi))

x11()
plot(res,asp=1,col=my.palette,pos.legend=3)
plot(db.regul.baliste)


########################
#moyenne pour biomasse 
abondance=NULL
biomasse=NULL
surface.presence=NULL
m=NULL
poids.moyen <-  485.05 #gramme
grille <- db.grid.init(db.regul.baliste,nodes=200)

for(i in 1:5){
variable <- db.extract(res,paste("Variable.",i,sep=""))
classe <- db.extract(res,"classe")

surface.presence[i] <- sum(classe==i,na.rm=T)*(grille@dx[1]*grille@dx[2])*60*60*1852^2

##moyenne de Sa quand incatrice ==1
m[i] = mean(db.regul.baliste.ind[,4][db.regul.baliste.ind[,9+i]==1])
abondance[i] <- m[i]*surface.presence/(1852^2*4*pi*10^(-39.3/10)) #en nbr d'individus
biomasse[i] <- abondance[i] * poids.moyen /10^6

}
sum(biomasse)
2825222*(grille@dx[1]*grille@dx[2])*60*60*1852^2


saveRDS(res, file = "res_kriMAF_baliste.rds")
readRDS(file = "res_kriMAF_baliste.rds")
res= res_kriMAF_baliste
