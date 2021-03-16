###
# Code for interpolating sA around Fernando de Noronha 500m Adri for conditional simulations
# Auteur : nicolas.bez@ird.fr & ju.salvetat@gmail.com  adrien.brunel.pro@protonmail.com
# Date : 16/03/2021
###


root.dir = getwd()

###################################################
#XXX
###################################################

  # loading of the functions used to unfold round dataset
  source(paste0(root.dir,"/2_functions/trapeze number 2021_03_03.txt"))
  source(paste0(root.dir,"/2_functions/transformation des coordonnées par quadrilatere 2021_03_03.txt"))
  source(paste0(root.dir,"/2_functions/projection point on line.txt"))



  # loading local functions (replication,extension,variogram computation)
  source(paste0(root.dir,"/2_functions/local_functions.R"))
  
  # loading local functions (replication,extension,variogram computation)
  my.polygone_shelf <- readRDS(paste0(root.dir,"/1_data/my.polygone_shelf.rds"))
  
  # Loading reference lines required for projecting datapoints.
  # References lines have been defined wrt to depth threshold equal to 100 m.
  # The two lines must have the same nb of points, but must not be convexe.
  # They must be provided as closed lines.
  # The quadrats must be trapezoides 
  source(paste0(root.dir,"/1_data/ref.line.R"))
  source(paste0(root.dir,"/1_data/ref.line.2.R"))


###################################################
# Read/process data
###################################################


  # Read --------------------------------------------------------------------
    my.dataframe <- read.csv(paste0(root.dir,"/1_data/F123_500m_Averaged_SaSum.csv"),sep = ",",header=TRUE)
  

  # Transform -----------------------------------------------------------------
    #excluding repeated magic square 
    my.dataframe$time_formatR = as.POSIXct((my.dataframe$Time_mean-719529)*86400,origin="1970-01-01",tz='UTC')
    carre_formatR <- strptime("2019-04-18 14:58:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
    carre_formatR_fin <- strptime("2019-04-20 10:00:00 UTC", "%Y-%m-%d %H:%M:%OS",tz='UTC')
    my.dataframe = my.dataframe[!(my.dataframe$time_formatR>carre_formatR & my.dataframe$time_formatR<carre_formatR_fin),]
    colnames(my.dataframe)
    # log
    my.dataframe$SaFish_logSa70 = log10(my.dataframe$SaFish_70+1)
    my.dataframe$SaFish_logSa200 = log10(my.dataframe$SaFish_200+1)
    # conservation of the four useful variables
    my.dataframe <- my.dataframe[,c(20,21,27,29)]
    names(my.dataframe) <- c('Latitude','Longitude','FAROFA','LogSa70')
    #excluding the star west of Fernando Drina submount
    sel <- my.dataframe$Longitude > -32.52
    my.dataframe <- my.dataframe[sel,]
    
    #hist(my.dataframe$SaFish_70)
    #hist(my.dataframe$SaFish_logSa70)

  
    #
    
    tmp <- my.polygone_shelf@sets[[1]]
    poly.Fernando <- polygon.create(tmp$x[c(1:25,48,49)],tmp$y[c(1:25,48,49)])
    tmp <- chull(my.dataframe$Longitude,my.dataframe$Latitude)
    poly.shelf <- list(x=my.dataframe$Longitude[tmp],y=my.dataframe$Latitude[tmp])
    poly.shelf <- polygon.create(x=poly.shelf$x,y=poly.shelf$y)



    
    
###################################################
# Read/process data
###################################################

  # DB for RGeostats
  db <- db.create(my.dataframe)
  db <- db.locate(db,c(3,2),loctype="x")
  db <- db.locate(db,c("LogSa70"),"z")
  hist(db$LogSa70)
  hist(db[,"LogSa70"],breaks=100,col="grey",main="",xlab="LogSa70")
  
  # display
  plot(db,col=my.palette,asp=1)
  plot(poly.Fernando,add=T)



###################################################
# Unfold/projection
###################################################
  # Projection  -------------------------------------------------------------
    tmp <- f.trapeze.nb(db$Longitude,db$Latitude,ref.line,ref.line.2)
    tmp2 <- proj.conformal(db$Longitude,db$Latitude,
                           ref.line = ref.line,
                           ref.line.2 = ref.line.2,tmp$trapeze.nb)
    tmp2$x1 <- tmp2$x1/20 # in order to get roughly similar scales along et across
    
    db.unfold <- db.create(tmp2)
    db.unfold <- db.locate(db.unfold,c(2,3),"x")
    db.unfold <- db.add(db.unfold,LogSa70=db$LogSa70)
    db.unfold <- db.add(db.unfold,tmp,loctype=NA)
    
    plot(db.unfold,col=my.palette)
    abline(v=0.05)
    # adding lines corresponding to the different trapeze
    yref=ref.line$y
    xref=ref.line$x
    nref=length(xref) 
    dy=yref[-1]-yref[-nref]
    dx=xref[-1]-xref[-nref]
    d = sqrt(dx^2+dy^2) 
    dd = c(0,cumsum(d)) 
    abline(h=dd)

  # Drift  ------------------------------------------------------------------
  x.breaks <- seq(0,0.1,0.001)
  x.mids <- (x.breaks[-1]+x.breaks[-length(x.breaks)])/2
  residuals <- db.unfold$LogSa70*NA
  
  x <- db.unfold$x1
  z <- db.unfold$LogSa70
  tmp.mean <- tapply(z,cut(x,breaks=x.breaks),FUN="mean")
  tmp.sd <- tapply(z,cut(x,breaks=x.breaks),FUN="sd")
  drift.m <- approxfun(x=x.mids,y=tmp.mean,rule=2)
  drift.sd <- approxfun(x=x.mids,y=tmp.sd,rule=2)
  
  residuals <- (z-drift.m(x))/drift.sd(x)
  
  db.unfold <- db.add(db.unfold,LogSa70.residuals = residuals)


# # Variography in the unfolded system -----------------------------------------------------
# # 
# 
# vario <- vario.calc(db.unfold,lag=0.0005,nlag=40)#,dirvect=c(0,90))
# x11()
# plot(vario,npairdw=T,flag.norm=F,add=T)
# model.vario <- model.auto(vario,struct=c(1,3,3,8),npairdw=T,wmode=2,
#                           auth.aniso = T,
#                           col=my.palette[c(1,16)],flag.noreduce = T)
  
  
  
  
###################################################
# Anamorphosis (Petitgas et al. 2017)
###################################################

  # Display the histogram
  hist(db.unfold[,"LogSa70"],breaks=100,col="grey",main="",xlab="LogSa70")
  
  # Define the anamorphosis model
  model.anam <- anam.fit(db.unfold,draw=TRUE)
  db.unfold <- anam.z2y(db.unfold,anam=model.anam)
  print(db.unfold,flag.stats=TRUE,names="Gaussian.LogSa70.residuals")
  db.unfold <- db.rename(db.unfold,name="Gaussian.LogSa70.residuals",newname="Yp")
  ycut <- round(qnorm(sum(db.extract(db.unfold,"LogSa70") == 0) / db.unfold$nech),5)
  Y <- db.extract(db.unfold,"Yp")
  Y[Y < ycut] <- ycut
  hist(Y)
  db.unfold <- db.replace(db.unfold,"Yp",Y)
  print(db.unfold ,flag.stats=TRUE,names="Yp")

  # Variogram of untracted variable (Yp)
  vario.Yp = vario.calc(db.unfold,lag=0.0005,nlag=40)#,dirvect=c(0,90))
  x11()
  plot(vario.Yp,npairdw=T,flag.norm=F,add=T)
  model.vario.Yp <- model.auto(vario.Yp,struct=c(1,3,3,8),npairdw=T,wmode=2,
                               auth.aniso = T,
                               col=my.palette[c(1,16)],flag.noreduce = T)
  
  # Variogram truncated variable (Y)
  n.H <- 40
  ycut <- round(qnorm(sum(db.extract(db.unfold,"LogSa70") == 0) / db$nech),5)
  vario.Y <- vario.trans.cut(vario.Yp,ycut,n.H)
  model.vario.Y <- model.auto(vario.Y,struc=melem.name(c(1,2,3)),draw=F)
  
  x11()
  plot(vario.Yp,npairdw=T,inches=0.05,col="black",ylim=c(0,1.2))
  plot(vario.Y,npairdw=T,inches=0.05,col="red",add=TRUE)
  plot(model.vario.Y,npairdw=T,add=T,col="red")
  plot(model.vario.Yp,npairdw=T,add=T,col="black")
  legend(x="bottomright",legend=c("Variogram of Yp","Variogram of Y"),
         lty=c(1,1),col=c("black","red"))

  # define interval limit
  Ymax <- db.extract(db.unfold,name="Yp",flag.compress=FALSE)
  Ymin <- db.extract(db.unfold,name="Yp",flag.compress=FALSE)
  Ymin[Ymin <= ycut] <- -10
  db.unfold<-db.add(db.unfold,Ymax)
  db.unfold<-db.locate(db.unfold,db.unfold$natt,"upper")
  db.unfold<-db.add(db.unfold,Ymin)
  db.unfold<-db.locate(db.unfold,db.unfold$natt,"lower")
  
  # A Gibbs sampler
  db.unfold <-gibbs(db = db.unfold, model = model.vario.Y, seed = 232132,
             nboot = 10, niter = 100, flag.norm=FALSE, percent=0,
             toleps = 1,
             radix = "Gibbs", modify.target = TRUE)
  
  db.unfold<-db.rename(db.unfold,"Gibbs.G1","Y")
  print(db.unfold,flag.stats=TRUE,names="Y")
  
  # Histograms
  hist(db.unfold[,"Yp"],breaks=100,xlim=c(-4,4),ylim=c(0,300),main="",xlab="Yp")
  hist(db.unfold[,"Y"],breaks=100,xlim=c(-4,4),ylim=c(0,300),main="",xlab="Y")
  
  # Variogram gibbs variable (unused later ...)
  vario.Yg  <- vario.calc(db.unfold,lag=0.0005,nlag=40)#,dirvect=c(0,90))
  x11()
  plot(vario.Yg,npairdw=TRUE,inches=0.05)
  plot(model.vario.Y,add=TRUE,col="red")

  
###################################################
# Projected grid
###################################################
  
  # definition de la grille de krigeage
  # la grille ne couvre que les donnees 
  # i.e. sans les extensions qui ne servent qu'? assurer la continuite spatiale
  # option 1 : grille definie en geographique puis projet?e  
  grid <- db.grid.init(db,nodes=200)
  grid@x0 <- grid@x0+runif(2,-0.00001,0.00001)
  grid <- db.polygon(grid,poly.shelf,combine = "and")
  grid <- db.polygon(grid,poly.Fernando,flag.out=T,combine = "and")
  # projection de la grille
  # extraction des points actifs de la grille
  tmp <-db.extract(grid,2:3,flag.compress = T)
  # projection ATTENTION l'objet resultat s'appelle grid.conform mais n'est plus une grille ? cause de la projection
  tmp2 <- f.trapeze.nb(tmp$x1,tmp$x2,ref.line,ref.line.2)
  tmp3 <- proj.conformal(tmp$x1,tmp$x2,ref.line,ref.line.2,tmp2$trapeze.nb)
  tmp3$x1 <- tmp3$x1/20
  db.pseudogrid.unfold <- db.create(tmp3)
  db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,tmp2,loctype=NA)
  db.unfold.extended <- f.extension.bilaterale(db.unfold)
  
  
###################################################
# Conditional simulation
###################################################

# Conditional simulation of the Gaussian variable
grid.simu <- simtub(dbin=db.unfold, dbout=db.pseudogrid.unfold, model=model.vario.Y,
                    neigh=neigh.create(type=0,ndim=2), uc = "", mean
                    = 0, seed = 232132,
                    nbsimu = 1, nbtuba = 1000,
                    radix = "Simu", modify.target = TRUE)
grid.simu <- db.rename(grid.simu,"Simu.Y.S1","Simu.Y")
print(grid.simu,flag.stats=TRUE,names="Simu.Y")
# Display
plot(grid.simu,name="Simu.Y",
     asp=1/cos(mean(db.extract(db=db.unfold,names="x2"))*pi/180),
     pos.legend=5,flag.proj=F,title="")
map("worldHires",add=T,fill=T,col=8)
# Transform gaussian conditional simulation into raw scale
grid.simu <- anam.y2z(grid.simu,name="Simu.Y",anam=model.anam)
# Display
plot(grid.simu,name="Raw.Simu.Y",
     asp=1/cos(mean(db.extract(db=db.unfold,names="x2"))*pi/180),
     pos.legend=5,flag.proj=F,title="",zlim=c(10,18000))
map("worldHires",add=T,fill=T,col=8)

# unfold 

zest <- rep(NA,db.pseudogrid.unfold$nech)
x <- db.pseudogrid.unfold$x1
z <- db.pseudogrid.unfold$Kriging.LogSa70.residuals.estim
zest <- z*drift.sd(x)+drift.m(x)
   
# db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,Kriging.Depth.estim = zest)

tmp <- rep(NA,db.pseudogrid.unfold$nech)
tmp[kri.shelf$sel] <- kri.shelf$Kriging.Depth.residuals.estim[kri.shelf$sel]
tmp[kri.shelfbreak$sel] <- kri.shelfbreak$Kriging.Depth.residuals.estim[kri.shelfbreak$sel]
kri.unfold <- db.pseudogrid.unfold
kri.unfold <- db.add(kri.unfold,Kriging.Depth.residuals.estim=tmp)



simu.unfold <- db.pseudogrid.unfold
simu.unfold <- db.add(simu.unfold,grid.simu[,8])
zest <- rep(NA,db.pseudogrid.unfold$nech)
x <- kri.unfold$x1
z <- kri.unfold$Kriging.LogSa70.residuals.estim
zest <- z*drift.sd(x)+drift.m(x)
#   
db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,Kriging.LogSa70.estim = zest)


tmp <- rep(NA,db.pseudogrid.unfold$nech)
tmp[kri.shelf$sel] <- kri.shelf$Kriging.Depth.residuals.estim[kri.shelf$sel]
tmp[kri.shelfbreak$sel] <- kri.shelfbreak$Kriging.Depth.residuals.estim[kri.shelfbreak$sel]
kri.unfold <- db.pseudogrid.unfold
kri.unfold <- db.add(kri.unfold,Kriging.Depth.residuals.estim=tmp)

# Adding back the standardized drift 
zest <- rep(NA,db.pseudogrid.unfold$nech)
x <- kri.unfold$x1
z <- kri.unfold$Kriging.Depth.residuals.estim
zest <- z*drift.sd(x)+drift.m(x)
db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,Kriging.Depth.estim = zest)

# R?-affectation des valeurs krigees aux noeuds de la grille de depart dans l'espace geographique pour avoir une representation correcte
kri <- grid
sel <- grid[,5]
tmp <- rep(NA,grid$nech)
tmp[sel] <- db.pseudogrid.unfold$Kriging.LogSa70.estim
#tmp [tmp < 0] <-0
kri <- db.add(kri, K.estim = tmp)

plot(kri,name.image="K.estim",col=my.palette,pos.legend=3,zlim=c(0,100),asp=1)
plot(poly.Fernando,add=T)
plot(db,add=T,col=1,inches=0.05)

plot(kri,name.contour="K.estim",col=my.palette,pos.legend=3,levels=c(seq(0,100,10),200,500,1000))#zlim=c(0,100))










# Neighborhood ------------------------------------------------------------


neigh <- neigh.create(type=2,radius=c(0.02,0.01),
                      flag.continuous = T, dist.cont=0.5,
                      #flag.sector = T,nsect = 6,nsmax=60,
                      nmini=100,nmaxi=100,
                      flag.aniso=T)

# neigh.shelf <- neigh.create(type=2,radius=c(0.01,0.03),
#                               flag.continuous = T, dist.cont=0.5,
#                               #flag.sector = T,nsect = 6,nsmax=60,
#                               nmini=50,nmaxi=100,
#                               flag.aniso=T)
# neigh.shelfbreak <- neigh.create(type=2,radius=c(0.01,0.03),
#                               flag.continuous = T, dist.cont=0.5,
#                               #flag.sector = T,nsect = 6,nsmax=60,
#                               nmini=80,nmaxi=100,
#                               flag.aniso=T)
# # neigh.conform.test <- neigh.test(regul.pt.conform.extended,db.pseudogrid.unfold,model=model.conform,neigh=neigh.conform)


# # test
# tmp2 <- f.trapeze.nb(x=-32.41346,y=-3.914363,ref.line)
# tmp3 <- proj.conformal(x=-32.41346,y=-3.914363,ref.line,ref.line.2,tmp2$trapeze.nb)
# tmp3$x1 <- tmp3$x1/20
# tmp3 <- db.create(tmp3)
# tmp3 <- db.add(tmp3,Depth.drift = drift(drift.est,tmp3$x1),loctype = "f1")
# 
# kri.test <- krigtest(regul.pt.conform.extended,tmp3,
#                      uc=c(1,"f1"),
#                      model=model.shelf,neigh=neigh.shelf)
# kri.wgt <- kri.test$wgt[1:200,1]
# sel <- kri.wgt > 0
# symbols(kri.test$xyz[,1][sel],kri.test$xyz[,2][sel],circles = kri.wgt[sel],asp=1,inches = 0.1,col=1)
# symbols(kri.test$xyz[,1][!sel],kri.test$xyz[,2][!sel],circles = abs(kri.wgt[!sel]),add=T,inches = 0.1,fg=2)
# 
# points(tmp3$x1,tmp3$x2,col=2,pch=20)
# sum(kri.wgt*kri.test$data[1:100,1])

# Kriging -----------------------------------------------------------------

# duplication of the data before and after to avoid border effects

db.unfold.extended <- f.extension.bilaterale(db.unfold)

# # krigeage de l'ensemble des donn?es 
 db.pseudogrid.unfold <- kriging(db.unfold.extended,db.pseudogrid.unfold, flag.std = F, model=model.vario,neigh=neigh)
 kri.unfold <- db.pseudogrid.unfold
 kri.unfold <- db.add(kri.unfold,Kriging.LogSa70.residuals.estim=tmp)
 zest <- rep(NA,db.pseudogrid.unfold$nech)
 x <- kri.unfold$x1
 z <- kri.unfold$Kriging.LogSa70.residuals.estim
 zest <- z*drift.sd(x)+drift.m(x)
#   
 db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,Kriging.LogSa70.estim = zest)


tmp <- rep(NA,db.pseudogrid.unfold$nech)
tmp[kri.shelf$sel] <- kri.shelf$Kriging.Depth.residuals.estim[kri.shelf$sel]
tmp[kri.shelfbreak$sel] <- kri.shelfbreak$Kriging.Depth.residuals.estim[kri.shelfbreak$sel]
kri.unfold <- db.pseudogrid.unfold
kri.unfold <- db.add(kri.unfold,Kriging.Depth.residuals.estim=tmp)

# Adding back the standardized drift 
zest <- rep(NA,db.pseudogrid.unfold$nech)
x <- kri.unfold$x1
z <- kri.unfold$Kriging.Depth.residuals.estim
zest <- z*drift.sd(x)+drift.m(x)
db.pseudogrid.unfold <- db.add(db.pseudogrid.unfold,Kriging.Depth.estim = zest)

# R?-affectation des valeurs krigees aux noeuds de la grille de depart dans l'espace geographique pour avoir une representation correcte
kri <- grid
sel <- grid[,5]
tmp <- rep(NA,grid$nech)
tmp[sel] <- db.pseudogrid.unfold$Kriging.LogSa70.estim
#tmp [tmp < 0] <-0
kri <- db.add(kri, K.estim = tmp)

plot(kri,name.image="K.estim",col=my.palette,pos.legend=3,zlim=c(0,100),asp=1)
plot(poly.Fernando,add=T)
plot(db,add=T,col=1,inches=0.05)

plot(kri,name.contour="K.estim",col=my.palette,pos.legend=3,levels=c(seq(0,100,10),200,500,1000))#zlim=c(0,100))



bathy = kri@items[["K.estim"]]
lon = kri@items[["x1"]]
lat = kri@items[["x2"]]


table_kri = cbind(lon,lat,bathy)
table_kri = as.data.frame(table_kri)
table_kri = na.omit(table_kri)

write.csv(table_kri,"table_kri_bathy.csv")

matkri= matrix(bathy,200,200)  
matlon= matrix(lon,200,200)  
matlat= matrix(lat,200,200) 
write.csv(matkri,"matkribathy.csv",row.names = F)
write.csv(matlon,"matlonbathy.csv",row.names = F)
write.csv(matlat,"matlatbathy.csv",row.names = F)