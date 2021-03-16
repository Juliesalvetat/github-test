############### zeros effect

# Pre-requisite
projec.toggle(0)
rg.load("Demo.herring.sa.scot.db.data","db.data")
rg.load("Demo.herring.sa.scot.poly.data","poly.data")
rg.load("Demo.herring.sa.scot.grid.simu","grid.simu")
projec.define(db=db.data)
# Display the histogram
hist(db.data[,"sa"],breaks=100,col="grey",main="",xlab="sa")
# Define the anamorphosis model (Right figure)
model.anam <- anam.fit(db.data,type="emp",ndisc=db.data$nech,
                       sigma2e=800,draw=TRUE,title="")
# Transform the data into Gaussian
db.data <- anam.z2y(db.data,anam=model.anam)
print(db.data,flag.stats=TRUE,names="Gaussian.sa")
db.data <- db.rename(db.data,name="Gaussian.sa",newname="Yp")
ycut <- round(qnorm(sum(db.extract(db.data,"sa") == 0) / db.data$nech),5)
Y <- db.extract(db.data,"Yp")
hist(Y)
Y[Y < ycut] <- ycut
db.data <- db.replace(db.data,"Yp",Y)
print(db.data,flag.stats=TRUE,names="Yp")
# Modeling Gaussian variable Y
n.H <- 50
vario.Yp <- vario.calc(db.data,lag=2.5,nlag=50)
vario.Y <- vario.trans.cut(vario.Yp,ycut,n.H)
model.vario.Y <- model.auto(vario.Y,struc=melem.name(c(1,2,3)),draw=F)
plot(vario.Yp,npairdw=T,inches=0.05,col="black",ylim=c(0,1.2))
plot(vario.Y,npairdw=T,inches=0.05,col="red",add=TRUE)
plot(model.vario.Y,add=T,col="red")
legend(x="bottomright",legend=c("Variogram of Yp","Variogram of Y"),
       lty=c(1,1),col=c("black","red"))
Ymax <- db.extract(db.data,name="Yp",flag.compress=FALSE)
Ymin <- db.extract(db.data,name="Yp",flag.compress=FALSE)
Ymin[Ymin <= ycut] <- -10
db.data<-db.add(db.data,Ymax)
db.data<-db.locate(db.data,db.data$natt,"upper")
db.data<-db.add(db.data,Ymin)
db.data<-db.locate(db.data,db.data$natt,"lower")
# A Gibbs sampler
db.data <-gibbs(db = db.data, model = model.vario.Y, seed = 232132,
                nboot = 10, niter = 100, flag.norm=FALSE, percent=0,
                toleps = 1,
                radix = "Gibbs", modify.target = TRUE)
db.data<-db.rename(db.data,"Gibbs.G1","Y")
print(db.data,flag.stats=TRUE,names="Y")
# Histograms
hist(db.data[,"Yp"],breaks=100,xlim=c(-4,4),ylim=c(0,300),main="",xlab="Yp")
hist(db.data[,"Y"],breaks=100,xlim=c(-4,4),ylim=c(0,300),main="",xlab="Y")
vario.Yg <- vario.calc(db.data,lag=2.5,nlag=50)
plot(vario.Yg,npairdw=TRUE,inches=0.05)
plot(model.vario.Y,add=TRUE,col="red")

