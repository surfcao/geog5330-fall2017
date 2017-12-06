library(adehabitatLT)

# The basic analysis element is 'ltraj' and it needs coordinates and the time stamp to construct a trajectory. The 'ltraj' precompute descriptive parameters of the trajectory including the time and distance lags between the relocations. 

data(whale)
plot(whale)
class(whale)

# The 'ltraj' can be constructed from a data frame

data(puechabonsp)
locs <- puechabonsp$relocs
locs <- as.data.frame(locs)
head(locs)

class(locs$Date)

da <- as.POSIXct(strptime(as.character(locs$Date),"%y%m%d"))

puech <- as.ltraj(xy = locs[,c("X","Y")], date = da, id = locs$Name)
puech

head(puech[[1]])
plot(puech)

# The trajeocty can be regular or irregular observed

is.regular(puech)

# plotltr can help investigate the observation intervals

plotltr(puech, "dt/3600/24")

# divide the trajectory into bursts if the interval is longer than a threshold value

foo <- function(dt) {
     return(dt> (100*3600*24))
 }

puech2 <- cutltraj(puech, "foo(dt)", nextr = TRUE)
puech2

burst(puech2)[3:4] <- c("Chou.1992", "Chou.1993")
puech2


# The function which.ltraj can also be used to identify the bursts satisfying a condition. For example, imagine that we want to identify the bursts where the distance between successive relocations was greater than 2000 metres at least once:

bu <- which.ltraj(puech2, "dist>2000")

#This data frame contains the ID, burst ID and relocation numbers satisfying the specified criterion. We can then extract the bursts satisfying this criterion:

(puech2[burst(puech2)%in%bu$burst])

plot(puech2[5])


# make the irregular observed trajectory into regular ones by filling in missing values

plotltr(puech2, "dt/3600/24")

refda <- strptime("00:00", "%H:%M")
refda

puech3 <- setNA(puech2, refda, 1, units = "day")
puech3 

plotltr(puech3)

# rounding the time to define a regular trajectory

data(ibexraw)

plotltr(ibexraw, "dt/3600")

refda <- strptime("2003-06-01 00:00", "%Y-%m-%d %H:%M", tz="Europe/Paris")
ib2 <- setNA(ibexraw, refda, 4, units = "hour")
ib2

plotltr(ib2, "dt/3600")

ib3 <- sett0(ib2, refda, 4, units = "hour")
ib3
plotltr(ib3, "dt/3600")


## alight the trajectory with the same duration

is.sd(ib3)

ib4 <- set.limits(ib3, begin = "2003-06-01 00:00", dur = 14, units = "day", pattern = "%Y-%m-%d %H:%M", tz="Europe/Paris")

ib4

is.sd(ib4)

di <- sd2df(ib4, "dist")

##########################
# Analysis of trajectory
##########################

## randomness of the mising values

runsNAltraj(ib4)

# This trajectory is regular. The bear was monitored during one month, with one relocation every 30 minutes. 

data(bear)
runsNAltraj(bear)
plotNAltraj(bear)

# conver the time stamped trajectory into non time-stampled ones
bearI <- typeII2typeI(bear)
plot(bearI)


## rediscretize this trajectory with constant step length of 500 metres:
bearIr <- redisltraj(bearI, 500)
bearIr

plot(bearIr)

# the geometrical properties can be studied by examining the distribution of the relative angles.

sliwinltr(bearIr, function(x) mean(cos(x$rel.angle)), type="locs", step=30)

cosrelangle <- cos(bearIr[[1]]$rel.angle)

head(cosrelangle)


## interpolate for by equal time intervals

data(porpoise)
porpoise2 <- redisltraj(na.omit(porpoise[1:3]), 86400, type="time")
plot(porpoise2[2])

#trajdyn(ib4)

## test autocorrelation of the parameters

wawotest(bear)

plotltr(bear, "dist")

sliwinltr(bear, function(x) mean(na.omit(x$dist)), 5*48, type="locs")

#ACF function for distance

acfdist.ltraj(bear, lag=5, which="dist")

# test the relative angels
testang.ltraj(bear, "relative")

acfang.ltraj(bear, lag=5)

# Segment of a spatial trajectory

data(porpoise)
gus <- porpoise[1]
gus
plot(gus)

acfdist.ltraj(gus, "dist", lag=20)
plotltr(gus, "dist")


(tested.means <- round(seq(0, 130000, length = 10), 0))
(limod <- as.list(paste("dnorm(dist, mean =", tested.means, ", sd = 5000)")))
mod <- modpartltraj(gus, limod)

mod

bestpartmod(mod)

(pm <- partmod.ltraj(gus, 4, mod))

plot(pm)

plotltr(gus, "dist")

tmp <- lapply(1:length(pm$ltraj), function(i){coul <- c("red","green","blue")[as.numeric(factor(pm$stats$mod))[i]]
    lines(pm$ltraj[[i]]$date, rep(tested.means[pm$stats$mod[i]], nrow(pm$ltraj[[i]])),col=coul, lwd=2)})


# Simulate a null model of trajectory

sim <- simm.crw(1:1000, r=0.95)
plot(sim, addp=F)

sim <- simm.levy(x=1:1000)
plot(sim, addp=F)


# Null model

data(puechcirc)
data(puechabonsp)
xo <- coordinates(puechabonsp$map)
boar1 <- puechcirc[1]
plot(boar1, spixdf=puechabonsp$map, xlim=range(xo[,1]), ylim=range(xo[,2]))

plotfun <- function(x, par) {
 image(par)
     points(x[,1:2], pch=16)
     lines(x[,1:2])
 return(x) }


confun <- function(x, par) {
     ## Define a SpatialPointsDataFrame from the trajectory
     coordinates(x) <- x[,1:2]
     ## overlap the relocations x to the elevation map par
     jo <- join(x, par)
     ## checks that there are no missing value
     res <- all(!is.na(jo))
     ## return this check
     return(res) }

 nmo2 <- NMs.randomShiftRotation(na.omit(boar1), rshift = TRUE, rrot = TRUE,
                                 rx = range(xo[,1]), ry = range(xo[,2]),
                                 treatment.func = plotfun,
                                 treatment.par = puechabonsp$map[,1],
                                 constraint.func = confun,
                                 constraint.par = puechabonsp$map[,1], nrep=9)

set.seed(90909)
par(mfrow=c(3,3), mar=c(0,0,0,0))
resu <- testNM(nmo2, count = FALSE)

# compare the variance of the elevation. 

varel <- function(x, par) {
     coordinates(x) <- x[,1:2]
     jo <- join(x, par) return(var(jo)) }

nmo3 <- NMs.randomShiftRotation(na.omit(boar1), rshift = TRUE, rrot = TRUE,
                                 rx = range(xo[,1]), ry = range(xo[,2]),
                                 treatment.func = varel,
                                 treatment.par = puechabonsp$map[,1],
                                 constraint.func = confun,
                                 constraint.par = puechabonsp$map[,1], nrep=99)

sim <- testNM(nmo3, count=FALSE)

(obs <- varel(na.omit(boar1)[[1]], puechabonsp$map[,1]))
(ran <- as.randtest(unlist(sim[[1]]), obs))
plot(ran)


## Null model with multiple animals

data(puechcirc) 
plot(puechcirc) 


(boar <- bindltraj(puechcirc))

nmo4 <- NMs.randomShiftRotation(na.omit(boar), rshift = TRUE, rrot = TRUE,
                                 rx = range(xo[,1]), ry = range(xo[,2]),
                                 treatment.func = varel,
                                 treatment.par = puechabonsp$map[,1],
                                 constraint.func = confun,
                                 constraint.par = puechabonsp$map[,1], nrep=99)
sim2 <- testNM(nmo4, count=FALSE)

(obs <- lapply(na.omit(boar), function(x) {varel(x, puechabonsp$map[,1])}))

lapply(1:2, function(i) {as.randtest(unlist(sim2[[i]]), obs[[i]])})
