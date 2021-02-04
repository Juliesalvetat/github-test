library("plyr") 
library(RGeostats)
library(maps)
library(mapdata)
library(wesanderson)
# Definition d'une palette
my.palette = wes_palette("Zissou1", 16, type = "continuous")

# loading of the function used to unfold round dataset
source("C:/Users/jsalveta/Documents/1_WORK/Cours/geostat/transformation_des_coordonnees_par_trapeze 2020_02_07.txt")

#setwd("C:/Users/nbez/Documents/stages these/THESES/Theses en cours/Salvetat Julie/kriggingBathyFAROFA")


#my.dataframe <- read.csv(file = 'FAROFA123_logSa_250m.csv')
#my.dataframe <- read.csv(file = 'FAROFA123_logSa_1ping.csv')
my.dataframe <- read.csv(file = 'Sa_Horizontal_F123_1ping_FISH.csv')

my.dataframe = na.omit(my.dataframe)
my.dataframe$time_formatR = as.POSIXct((my.dataframe$Time-719529)*86400,origin="1970-01-01",tz='UTC')
carre_formatR <- strptime("2019-04-18 14:58:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
carre_formatR_fin <- strptime("2019-04-20 10:00:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
my.dataframe = my.dataframe[!(my.dataframe$time_formatR>carre_formatR & my.dataframe$time_formatR<carre_formatR_fin),]

my.dataframe$Horizontal_logSa70 = log10(my.dataframe$Horizontal_Sa70+1)
my.dataframe$Horizontal_logSa200 = log10(my.dataframe$Horizontal_Sa200+1)


# ATTENTION : Sa70 = Sa200 

# lecture du polygone 
my.polygone_shelf <- readRDS("my.polygone_shelf.rds")

# creation de 2 polygones : 1 pour l'ile 1, pour le shelfbreak.
tmp <- my.polygone_shelf@sets[[1]]
poly.Fernando <- polygon.create(tmp$x[c(1:25,48,49)],tmp$y[c(1:25,48,49)])
poly.shelf <- polygon.create(tmp$x[26:47],tmp$y[26:47])

# reduction de la resolution de Fernando pour creer un perimetre convexe clockwise ferm?
# la routine de unfolding se base en effet sur une ligne de r?ference convexe 
tmp <- chull(poly.Fernando@sets[[1]]$x[c(1:24)],poly.Fernando@sets[[1]]$y[c(1:24)])
sel.ref <- rev(tmp[order(tmp)]) # clockwise
sel.ref <- c(sel.ref,sel.ref[1]) # fermeture
ref.line <- list(x=poly.Fernando@sets[[1]]$x[sel.ref],y=poly.Fernando@sets[[1]]$y[sel.ref])

# debut de la ref.line en dehors de la zone dense pour recentrer graphiquement cette zone dense (?tape cosm?tique ...)
tmpx <- ref.line$x 
tmpy <- ref.line$y
n <- length(tmpx)
#tmpx <- tmpx[-n] ; tmpy <- tmpy[-n]
tmpxx <- c(tmpx[7:(n-1)],tmpx[1:6])
tmpyy <- c(tmpy[7:(n-1)],tmpy[1:6])
tmpxx[n] <- tmpxx[1] ;tmpyy[n] <- tmpyy[1]
ref.line <- list(x=tmpxx,y=tmpyy)

# creation DB
db <- db.create(my.dataframe)
db <- db.locate(db,c(2,3),loctype="x")
db <- db.locate(db,c("Horizontal_logSa70"),"z")
# Simplification des noms des variables
db <- db.rename(db,c("Horizontal_Sa70","Horizontal_logSa70","Bottom_Depth"),c("Sa70","LogSa70","Depth"))
sum(db$LogSa70==0)/db$nech
# ==> 67.8% de 0 !!!

# Regularisation
# Une fois acquis le fait que les differentes FAROFA peuvent etre melangees (resultats not shown), on fait des moyennes de Sa par pixel
# en melangeant les 3 FAROFA
regul.grid <- db.grid.init(db,nodes=400)
regul.grid <- db.stat.grid(db,regul.grid,names=c("LogSa70","Sa70","Depth","FAROFA"),radix="",modify.target = F)
regul.grid <- db.locate(regul.grid,"LogSa70","z")
regul.grid <- db.rename(regul.grid,c("x1","x2"),c("Longitude","Latitude"))

# caract?ristiques stat des donn?es r?gularisees
# Nbr de cellules actives r?duit ? 16 482 
sum(!is.na(regul.grid$LogSa70)) 
# atome en 0 r?duit ? 37%
sum(regul.grid$LogSa70==0,na.rm=T)/sum(!is.na(regul.grid$LogSa70)) 
# mais histogramme similaire though with lowest CV
hist(regul.grid$LogSa70, main="Histogram LogSa70 in regularised grid", xlab = "LogSa70")
hist(db$LogSa70)
# les cellules sont fortement mono farofa
hist(regul.grid$FAROFA)
sum(regul.grid$FAROFA==1 | regul.grid$FAROFA==2 |regul.grid$FAROFA==3,na.rm=T)/sum(!is.na(regul.grid$LogSa70))
##==> les cellules sont mono farofa ? 90%

# Ecriture des donnees regularisees sous forme de point et non de grille.
# Cela ppermet de faire des calculs de variogrammes avec tolerance sur la direction.
# vario.grid() ne calcule que le long des axes principaux de la grille et avec tous les points de la grille
# y compris toutes les valeurs NA !!!
regul.pt <- db.grid2point.sampling(regul.grid,npack=1,names=c("LogSa70","Sa70","Depth","FAROFA"))
regul.pt <- db.reduce(db.sel(regul.pt,!is.na(LogSa70)))
regul.pt <- db.locate(regul.pt,c("Depth","FAROFA"),NA)
plot(regul.pt,col=my.palette)

# projection des points 
tmp <- proj.conformal(regul.pt$Longitude,regul.pt$Latitude,ref.line = ref.line,orthogonal = F,closed=T)
trapeze.nb = res$trapeze.nb #ajout julie sauvegarde trapeze.nb avec ajout dans proj.conformal de assign("res",res,.GlobalEnv)
regul.pt.conform <- db.create(tmp)
regul.pt.conform  <- db.locate(regul.pt.conform,c(2,3),loctype="x") #ajout julie  sinon bug pas de dimension spatiale
regul.pt.conform <- db.add(regul.pt.conform,LogSa70=regul.pt[,4])
plot(regul.pt.conform,col=my.palette,asp=1) # l'argument asp=1 permet d'avoir des graphes "carr?s" i.e. distances repr?sent?es en x et en y similaires
abline(v=0)


#Analyse variographique de Log10Sa70 dans l'espace conforme ?tendu unilat?ralement
# dir = 0 : perpendicular to the coast
# dir = 90 : parallel to the coast

f.replication.unilaterale <- function(db,line=ref.line){
  yref=line$y ; xref=line$x
  nref=length(xref) # nbr de points dans la ligne de r?f?rence
  dd = c(0,cumsum(sqrt((xref[-1]-xref[-nref])^2+(yref[-1]-yref[-nref])^2))) ### longeur cumul?e
  # ajout d'une variable pour tracer les donnees ajoutees
  db <- db.add(db,code=1,auto.locate=T)
  # R?plication des donn?es ? la fin du d?ploiement
  tmp <- db.extract(db,names=db$names[-1])
  tmp$y <- tmp$y + dd[nref] #max(tmp$y[trapeze.fin]) #changemant de tmp$x2 a tmp$y
  tmp$code <- 2
  dbout <- db.append(db,tmp)
  #
  dbout
}

# variogramme des donn?es depliees
vario <- vario.calc(regul.pt.conform,lag=c(0.0004,0.0012),nlag=c(150,50),dirvect = c(0,90)) 

# variogramme // ? la cote entre les donnees depliees et leur replication
# en prenant soin que la variable codante soit differente pour eviter les doublons avec le calcul precedent
# puis ajout de ce variogramme compl?mentaire au pr?c?dent 
tmp <- f.replication.unilaterale(regul.pt.conform)
#
vario.tmp <- vario.calc(tmp,lag=0.0012,nlag=50,dirvect = 90,tolang = 45,tolcode=2,opt.code=2) 
#
v1 <- vario ; v2 <- vario.tmp ; v <- vario 
v@vardirs[[2]]$gg <- (v1@vardirs[[2]]$gg*v1@vardirs[[2]]$sw +v2@vardirs[[1]]$gg*v2@vardirs[[1]]$sw)/(v1@vardirs[[2]]$sw +v2@vardirs[[1]]$sw)
v@vardirs[[2]]$sw <- v1@vardirs[[2]]$sw +v2@vardirs[[1]]$sw 

plot(vario,npairdw=T,col=my.palette[c(1,16)])
plot(v,npairdw=T,col=my.palette[c(1,16)],add=T)
# tr?s peu de modification du variogramme initial mais le calcul est plus propre ...

vario <- v


# reduction du vario aux distances < 0.03 avant modelisation pour que le mod?le s'ajuste au mieux sur ces distances
# malheureusement il n'existe pas de fonction toute pr?te. il faut le faire ? la main.
# ou refaire le calcul de variogramme avec moins de pas de calcul (mais ca prend du temps ...)
tmp <- vario
for(i in 1:2){
sel <- tmp@vardirs[[i]]$hh < 0.03
tmp@vardirs[[i]]$hh <- tmp@vardirs[[i]]$hh[sel]
tmp@vardirs[[i]]$gg <- tmp@vardirs[[i]]$gg[sel]
tmp@vardirs[[i]]$sw <- tmp@vardirs[[i]]$sw[sel]
tmp@vardirs[[i]]@size <- sum(sel)
tmp@vardirs[[i]]@npas <- sum(sel)
tmp@vardirs[[i]]@npatot <- sum(sel)
}
vario.0.03 <- tmp

plot(vario.0.03,npairdw=T,col=my.palette[c(1,16)])
#
model.conform <- model.auto(vario.0.03,struct=c(1,2,2,2),npairdw=T,wmode=2,col=my.palette[c(1,16)])

### krgigeage

# creation d'un jeu de donn?es ?tendu afin d'?viter les effets de bord en cokrigeage
# on recopie en haut en bas les n.replicat trapezes adjacents

f.extension.bilaterale <- function(db,line=ref.line,n.replicat=4){
  # n.replicat : nbr de trapezes replicat en haut et en bas
  yref=line$y ; xref=line$x
  nref=length(xref) # nbr de points dans la ligne de r?f?rence
  dd = c(0,cumsum(sqrt((xref[-1]-xref[-nref])^2+(yref[-1]-yref[-nref])^2))) ### longeur cumul?e
  # ajout d'une variable pour tracer les donnees ajoutees
  db <- db.add(db,original=T,type.locate=F)
  # Ajout des trapezes ? la fin du d?ploiement
  tmp <- db.extract(db.reduce(db.sel(db,trapeze.nb <=n.replicat)),names=db$names[-1])
  tmp$y <- tmp$y + dd[nref] #max(tmp$y[trapeze.fin])
  tmp$original <- F
  dbout <- db.append(db,tmp)
  # Ajout des trapezes avant le d?but du d?ploiement en inversant la coordonn?es y
  tmp <- db.extract(db.reduce(db.sel(db,trapeze.nb >=(max(trapeze.nb)-(n.replicat-1)))),names=db$names[-1])
  tmp$y <- tmp$y - dd[nref] #max(tmp$y[trapeze.fin]) # pareil changement tmp$x2 a tmp$y
  tmp$original <- F
  dbout <- db.append(dbout,tmp)
  #
  dbout
}

regul.pt.conform.extended <- f.extension.bilaterale (regul.pt.conform)


# definition de la grille de krigeage
# la grille ne couvre que les donnees 
# i.e. sans les extensions qui ne servent qu'? assurer la continuite spatiale
# option 1 : grille definie en geographique puis projet?e  
grid <- db.grid.init(regul.pt,nodes=150)
grid@x0 <- grid@x0+runif(2,-0.001,0.001)
grid <- db.polygon(grid,poly.shelf,combine = "and")
grid <- db.polygon(grid,poly.Fernando,flag.out=T,combine = "and")
# projection de la grille
# extraction des points actifs de la grille
tmp <-db.extract(grid,2:3,flag.compress = T)
# projection ATTENTION l'objet resultat s'appelle grid.conform mais n'est plus une grille ? cause de la projection
tmp <- proj.conformal(tmp$x1,tmp$x2,ref.line,orthogonal = F,closed=T)
#trapeze.nb = res$trapeze.nb
db.pseudogrid.conform <- db.create(tmp)
db.pseudogrid.conform <- db.locate(db.pseudogrid.conform,c(2,3),loctype="x")

plot(db.pseudogrid.conform)

#voisinage
#neigh <- neigh.create(ndim=2,nmini=20,nmaxi = 100,radius=0.02)
rayon <- 0.02
neigh.conform <- neigh.create(type=2,radius=c(0.02,0.02),flag.continuous = T,
                              flag.sector = T,nsect = 6,nsmax=40,nmini=10,nmaxi=240,
                              flag.aniso=T)
#neigh.conform.test <- neigh.test(regul.pt.conform.extended,db.pseudogrid.conform,model=model.conform,neigh=neigh.conform)

# attention avec 150 noeuds dans la definition de la grille ca prend 5 minutes 
# on krige vers la grille cal?e sur les donn?es mais avec des donn?es ?tendues

# kri.conform <- kriging(f.extension.bilaterale(regul.pt.conform),db.pseudogrid.conform,model=model.conform,neigh=neigh.conform,ndisc = 50)#,calcul="block")
# kri.conform <- kriging(f.extension.bilaterale(regul.pt.conform),grid,model=model.conform,neigh=neigh.conform,ndisc = 50,calcul="block")
# 
# db.pseudogrid.conform


#grid.kri <-  db.grid.init(f.extension.bilaterale(regul.pt.conform),nodes=150)
#grid.kri  <- db.locate(regul.grid,"LogSa70","z")
#db.grid.conform <- db.grid.init(db.pseudogrid.conform)
#kri.conform <- kriging(f.extension.bilaterale(regul.pt.conform),db.grid.conform,model=model.conform,neigh=neigh.conform,calcul="block")
#x11()
#plot(kri.conform,name.image="Kriging.LogSa70.estim",col=my.palette,pos.legend=3)


#kri.conform <- kriging(f.extension.bilaterale(regul.pt.conform),db.pseudogrid.conform,model=model.conform,neigh=neigh.conform,calcul="block")


# R?-affectation des valeurs krigees aux noeuds de la grille de depart dans l'espace geographique pour avoir une representation correcte
kri <- grid
sel <- grid[,5]
tmp <- rep(NA,grid$nech)
tmp[sel] <- kri.conform$Kriging.LogSa70.estim
tmp [tmp < 0] <-0
kri <- db.add(kri, K.estim = tmp)
tmp[sel] <- kri.conform$Kriging.LogSa70.stdev
kri <- db.add(kri,K.stdev = tmp)
x11()
plot(kri,name.image="K.estim",col=my.palette,pos.legend=3)
plot(poly.Fernando,add=T)

#######################
#######################
# Passage aux indicatrices
#######################
#######################

# Ajout des indicatrices ? la base de donn?es regularisees
# Ind1 = presence/absence
# Ind2-4 = quartile des donn?es positives
dbin <- regul.pt
zcut <- quantile(dbin$LogSa70[dbin$LogSa70 > 0],seq(0.25,1,0.25),na.rm=T)
#my.limits <- limits.create(mini=c(0,zcut[-4]),maxi=zcut,flag.zcut.int=flag.zcut.int,incmini=T,incmaxi=c(F,F,F,T))
my.limits <- limits.create(mini=c(0,0.000001,zcut[-4]),maxi=c(0.000001,zcut),incmini=c(T,rep(F,4)),incmaxi=rep(T,5),flag.zcut.zero=T)
print(my.limits)
dbout <- db.indicator(dbin, limits=my.limits,radix="Ind")
regul.pt.ind <-  dbout

# projection des points 
tmp <- proj.conformal(regul.pt.ind$Longitude,regul.pt.ind$Latitude,ref.line = ref.line,orthogonal = F,closed=T)
regul.pt.ind.conform <- db.create(tmp)
regul.pt.ind.conform <- db.add(regul.pt.ind.conform,regul.pt.ind[,c(4,8:12)],loctype = NA)
regul.pt.ind.conform <- db.locate(regul.pt.ind.conform,c(2,3),loctype="x")

# representation geographique des indicatrices de classe
plot(db.sel(regul.pt.ind,Ind.LogSa70.1==1),inches=0.1,pch=20,col=my.palette[1])
plot(db.sel(regul.pt.ind,Ind.LogSa70.2==1),inches=0.1,pch=20,col=my.palette[5],add=T)
plot(db.sel(regul.pt.ind,Ind.LogSa70.3==1),inches=0.1,pch=20,col=my.palette[9],add=T)
plot(db.sel(regul.pt.ind,Ind.LogSa70.4==1),inches=0.1,pch=20,col=my.palette[13],add=T)
plot(db.sel(regul.pt.ind,Ind.LogSa70.5==1),inches=0.1,pch=20,col=my.palette[16],add=T)

plot(poly.Fernando,add=T)

regul.pt.ind.conform <- db.locate(regul.pt.ind.conform,names="Ind*","z")
# covariance des indicatrices
covnc.ind <- vario.calc(regul.pt.ind.conform,calcul="covnc",lag=c(0.0004,0.0012),nlag=c(70,25),dirvect = c(0,90))
plot(covnc.ind,col=my.palette[c(1,16)],npairdw=T)
# on voit des dissym?tries dans les variogrammes crois?s assez int?ressantes 
# elles racontent que la probabilit? de rencontrer du riche quand on est dans du moyen est plus forte en allant vers la cote que vers le large. 
# on pourra en reparler mais la mod?lisation par indicatrice permet de d?crire le mode d'imbrication spatiale des niveaux de densit?.

### MAF
# on travaille avec 3 indicatrices ; la 4eme se d?duisant des 3 premi?res.
#regul.pt.ind.conform <- db.locate(regul.pt.ind.conform,names=19:13, loctype = "z", flag.locnew = TRUE)
maf <- maf.calc(regul.pt.ind.conform,h0=0.0008,dh=0.0004)
regul.pt.ind.conform.maf <- pca.z2f(regul.pt.ind.conform,pca=maf,verbose = F)

# variogramme des MAFs
vario.maf <- vario.calc(regul.pt.ind.conform.maf,lag=c(0.0004,0.0012),nlag=c(70,25),dirvect = c(0,90)) 
plot(vario.maf,npairdw=T,col=my.palette[c(1,16)])

# extraction des variogrammes simples en vue de faire des Krigeages Simples
for(i in 1:5) assign(paste("vario.maf.",i,sep=""),vario.reduce(vario.maf,varcols=i,flag.vario=T))
for(i in 1:5) plot(get(paste("vario.maf.",i,sep="")),add=i>1,ylim=c(0,1.4),npairdw=T,col=my.palette[c(1,16)])

for(i in 1:5){
  for(j in (i+1):5){
    plot(vario.maf,varcols=i,varcols2=j,add=!(i==1 & j==2),ylim=c(-0.1,0.1),xlim=c(0,0.01),npairdw=T,col=my.palette[c(1,16)])
  }
}

# reduction du vario du MAF3 aux distances < 0.0005 avant modelisation pour que le mod?le s'ajuste au mieux sur ces distances
tmp <- vario.maf.3
for(i in 1:2){
  sel <- tmp@vardirs[[i]]$hh < 0.005
  tmp@vardirs[[i]]$hh <- tmp@vardirs[[i]]$hh[sel]
  tmp@vardirs[[i]]$gg <- tmp@vardirs[[i]]$gg[sel]
  tmp@vardirs[[i]]$sw <- tmp@vardirs[[i]]$sw[sel]
  tmp@vardirs[[i]]@size <- sum(sel)
  tmp@vardirs[[i]]@npas <- sum(sel)
  tmp@vardirs[[i]]@npatot <- sum(sel)
}
vario.maf.3 <- tmp

# modelisation des variogrammes
for(i in 1:5) assign(paste("model.maf.",i,sep=""),model.auto(get(paste("vario.maf.",i,sep="")),struct=c(1,2,2,2),wmode=2,npairdw=T))

# r?alisation des krigeages monovariables des 3 MAFs*
# adaptation de la taille du voisinage ? la port?e des mod?les
# attention si le voisinage est trop restreint il apparait des zones non estim?es ...
rayon <- c(0.02,0.01,0.01,0.01,0.01)

# temps de clacul des 3 krigeages monovariable environ 10-20 mn.
# en cokrigeage multivari? la nuit ne suffit pas ...
for(i in 1:5){
  cat("kriging MAF ",i,"\n")
  neigh.conform <- neigh.create(type=2,radius=rayon[i],flag.continuous = T,
                                flag.sector = T,nsect = 6,nsmax=40,nmini=10,nmaxi=240,
                                flag.aniso=T)
  regul.pt.ind.conform.maf <- db.locate(regul.pt.ind.conform.maf,paste("Factor.",i,sep=""),"z", flag.locnew = TRUE)
  res <- kriging(f.extension.bilaterale(regul.pt.ind.conform.maf),db.pseudogrid.conform,model=get(paste("model.maf.",i,sep="")),neigh=neigh.conform)
  assign(paste("kri.maf.",i,sep=""),res)
}

# R?-affectation des valeurs krigees aux noeuds de la grille de depart dans l'espace geographique pour avoir une representation correcte
kri <- grid
sel <- grid[,5]
tmp <- rep(NA,grid$nech)
for(i in 1:5){
  tmp0 <- get(paste("kri.maf.",i,sep=""))
  tmp[sel] <- tmp0[,db.getcols(tmp0,"z",1)]
  kri <- db.add(kri, tmp)
  kri <- db.rename(kri,"tmp",paste("K.estim.maf.",i,sep=""))
}

plot(kri,name.image="*maf.1",col=my.palette,pos.legend=3)
plot(poly.Fernando,add=T)

### Reconstruction des densit?s interpol?es ? partir des MAFs interpol?s
kri <- db.locate(kri,6:10,"z")
res <- pca.f2z(kri,pca=maf,verbose = F)

# seuillage des krigeages d'indicatrices ? [0,1]
for(i in 11:15){
  res[,i][res[,i] > 1] <- 1
  res[,i][res[,i] < 0] <- 0
}

plot(res,name=12,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(LogSa < q25+)")
plot(res,name=13,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(q25+ < LogSa < q50+)")
plot(res,name=14,col=my.palette,pos.legend=3,zlim=c(0,1),title="P(q50+ < LogSa < q75+)",asp=1,las=1) ; plot(poly.Fernando,add=T)
plot(res,name=15,col=my.palette,pos.legend=3,zlim=c(0,1),title="P( LogSa > q75+)",asp=1,las=1) ; plot(poly.Fernando,add=T)

res <- db.compare(res,fun="maxi",names=11:15)
res <- db.add(res, classe = 1*(Variable.1==maxi)+2*(Variable.2==maxi)+3*(Variable.3==maxi)+4*(Variable.4==maxi)+5*(Variable.5==maxi))

x11()
plot(res_kriMAF,asp=1,col=my.palette,pos.legend=3)
plot(db.regul.baliste)

saveRDS(res, file = "res_kriMAF.rds")
readRDS(file = "res_kriMAF.rds")
res =res_kriMAF

surface.presence=NULL
m=NULL
for(i in 1:5){
  variable <- db.extract(res,paste("Variable.",i,sep=""))
  classe <- db.extract(res,"classe")
  surface.presence[i] <- sum(classe==i,na.rm=T)*(grid@dx[1]*grid@dx[2])*60*60*1852^2
  m[i] = mean(regul.pt.ind[,5][regul.pt.ind[,8+i]==1])
}

table_surface_mean = cbind(surface.presence,m)
write.csv(table_surface_mean,"table_surface_mean_kriMAF.csv")


varkri = res@items[["Kriging"]]
lon = res@items[["x1"]]
lat = res@items[["x2"]]
classe = res@items[["classe"]]

# matkri= matrix(varkri,150,150)  
# matlon= matrix(lon,150,150)  
# matlat= matrix(lat,150,150)  

table_kri = cbind(lon,lat,classe)
table_kri = as.data.frame(table_kri)
table_kri = na.omit(table_kri)

write.csv(table_kri,"table_kriMAF.csv")
