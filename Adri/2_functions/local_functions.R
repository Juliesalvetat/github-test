# Definitions of ad hoc functions --------------------------------------------------------

# Replication of the alined database 
# Used to compute clocwise and anticlockwise distances
f.replication.unilaterale <- function(db,line=ref.line){
  yref=line$y ; xref=line$x
  nref=length(xref) # nbr de points dans la ligne de r?f?rence
  dd = c(0,cumsum(sqrt((xref[-1]-xref[-nref])^2+(yref[-1]-yref[-nref])^2))) ### longeur cumul?e
  # ajout d'une variable pour tracer les donnees ajoutees
  db <- db.add(db,code=1,auto.locate=T)
  # R?plication des donn?es ? la fin du d?ploiement
  tmp <- db.extract(db,names=db$names[-1])
  tmp$x2 <- tmp$x2 + dd[nref] #max(tmp$y[trapeze.fin])
  tmp$code <- 2
  dbout <- db.append(db,tmp)
  #
  dbout
}


# Duplication of the alined database on both sides to avoid border effects when kriging
f.extension.bilaterale <- function(db,line=ref.line,n.replicat=4){
  # n.replicat : nbr de trapezes replicat en haut et en bas
  yref=line$y ; xref=line$x
  nref=length(xref) # nbr de points dans la ligne de r?f?rence
  dd = c(0,cumsum(sqrt((xref[-1]-xref[-nref])^2+(yref[-1]-yref[-nref])^2))) ### longeur cumul?e
  # ajout d'une variable pour tracer les donnees ajoutees
  db <- db.add(db,original=T,type.locate=F)
  # Ajout des trapezes ? la fin du d?ploiement
  tmp <- db.extract(db.reduce(db.sel(db,trapeze.nb <=n.replicat)),names=db$names[-1])
  tmp$x2 <- tmp$x2 + dd[nref] #max(tmp$y[trapeze.fin])
  tmp$original <- F
  dbout <- db.append(db,tmp)
  # Ajout des trapezes avant le d?but du d?ploiement en inversant la coordonn?es y
  tmp <- db.extract(db.reduce(db.sel(db,trapeze.nb >=(max(trapeze.nb)-(n.replicat-1)))),names=db$names[-1])
  tmp$x2 <- tmp$x2 - dd[nref] #max(tmp$y[trapeze.fin])
  tmp$original <- F
  dbout <- db.append(dbout,tmp)
  #
  dbout
}

### calcul rigoureux du variogrammes avec calcul rigoureux des distances // et des
### distance orthogonale a la cote
f.vario.calc <- function(mydb,lag=1, nlag=10){
  # db are supposed to have x1 and x2 coordinates names and to have a variable with locator z1
  # if my db gets an active selection only active points must be provided (mydb$x1 takes all the datapoints) 
  breaks <- seq(lag/2, by = lag,length=nlag+1)
  # distance orthogonale
  # et suppression du triangle inf?rieur pour all?ger les calculs finaux
  D.ortho <- abs(outer(mydb$x1,mydb$x1,FUN="-"))
  D.ortho <- D.ortho[upper.tri(D.ortho)]
  # Distance parallele 
  # Distance clockwise
  D.par.cw <- abs(outer(mydb$x2,mydb$x2,FUN="-"))
  # Distance parallele anti-clockwise
  tmp <- f.replication.unilaterale(mydb)
  tmp.x2 <- tmp$x2
  tmp.code <- tmp$code
  D.par.acw <- abs(outer(tmp.x2[tmp.code==1],tmp.x2[tmp.code==2],FUN="-"))
  # Distance parallele
  D.par <- pmin(D.par.cw,D.par.acw)
  D.par <- D.par[upper.tri(D.par)]
  # distance euclidienne
  D <- sqrt(D.par^2+D.ortho^2)
  # angle et selection des directions de calcul
  alpha <- atan(D.par/D.ortho)*180/pi
  dir.ortho <- alpha < 45
  dir.par <- alpha >= 45
  # 1/2 square dif of depth
  tmp <- mydb[,db.getcols(mydb,"z",1)]
  sigma2 <- var(tmp)
  D.z2 <- 0.5*outer(tmp,tmp,FUN="-")^2
  D.z2 <- D.z2[upper.tri(D.z2)]
  # variogramme othogonal
  tmp <- cut(D,breaks)
  g.ortho <- aggregate(D.z2[dir.ortho]~tmp[dir.ortho],FUN="mean")[,2]
  dist.ortho <- aggregate(D[dir.ortho]~tmp[dir.ortho],FUN="mean")[,2]
  npairs.ortho <- aggregate(D[dir.ortho]~tmp[dir.ortho],FUN="length")[,2]
  g.par <- aggregate(D.z2[dir.par]~tmp[dir.par],FUN="mean")[,2]
  dist.par <- aggregate(D[dir.par]~tmp[dir.par],FUN="mean")[,2]
  npairs.par <- aggregate(D[dir.par]~tmp[dir.par],FUN="length")[,2]
  # construction du resultat
  res <- vario.create(lag=lag,gg=g.ortho,hh=dist.ortho,sw=npairs.ortho,toldis=0.5,tolang=45,vars=sigma2)
  res <- vario.create(lag=lag,gg=g.par,hh=dist.par,sw=npairs.par,toldis=0.5,tolang=45,vario=res,codir=c(0,1))
  res
}