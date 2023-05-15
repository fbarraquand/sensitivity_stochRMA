### R file for plotting power spectrum

RMA = read.table("RMA.txt")
ncol(RMA)
names(RMA) = c("n","p","n2","p2","time","K")

RMA_I = read.table("RMA_I.txt")
ncol(RMA_I)
names(RMA_I) = c("n","p","n2","p2","time","K")

RMA_T = read.table("RMA_T.txt")
ncol(RMA_T)
names(RMA_T) = c("n","p","n2","p2","time","K")

## Selecting a couple of K close to 0.1, 0.5, 1.0 
RMA[RMA$K %in% c(0.128,0.520,1.024),]

plot(RMA$time[RMA$K==0.128],RMA$n[RMA$K==0.128],type="o",xlab = "Time",ylab = "Prey abundance")

par(mfrow=c(3,3))
spectrum(RMA$n[RMA$K==0.128])
spectrum(RMA_I$n[RMA_I$K==0.128])
spectrum(RMA_T$n[RMA_T$K==0.128])

spectrum(RMA$n[RMA$K==0.520])
spectrum(RMA_I$n[RMA_I$K==0.520])
spectrum(RMA_T$n[RMA_T$K==0.520])

spectrum(RMA$n[RMA$K==1.024])
spectrum(RMA_I$n[RMA_I$K==1.024])
spectrum(RMA_T$n[RMA_T$K==1.024])

### ACF -- do we see some periodicity here? 
par(mfrow=c(3,3))
acf(RMA$n[RMA$K==0.128],lag.max = 250)
acf(RMA_I$n[RMA_I$K==0.128],lag.max = 250)
acf(RMA_T$n[RMA_T$K==0.128],lag.max = 250)

acf(RMA$n[RMA$K==0.520],lag.max = 250)
acf(RMA_I$n[RMA_I$K==0.520],lag.max = 250)
acf(RMA_T$n[RMA_T$K==0.520],lag.max = 250)

acf(RMA$n[RMA$K==1.024],lag.max = 250)
acf(RMA_I$n[RMA_I$K==1.024],lag.max = 250)
acf(RMA_T$n[RMA_T$K==1.024],lag.max = 250)
## Seems like some long-run BUT these can perhaps arise out of randomness?

## For 50 timesteps -- given we cannot test for much more
par(mfrow=c(3,3))
acf(RMA$n[RMA$K==0.128],lag.max = 50)
acf(RMA_I$n[RMA_I$K==0.128],lag.max = 50)
acf(RMA_T$n[RMA_T$K==0.128],lag.max = 50)

acf(RMA$n[RMA$K==0.520],lag.max = 50)
acf(RMA_I$n[RMA_I$K==0.520],lag.max = 50)
acf(RMA_T$n[RMA_T$K==0.520],lag.max = 50)

acf(RMA$n[RMA$K==1.024],lag.max = 50)
acf(RMA_I$n[RMA_I$K==1.024],lag.max = 50)
acf(RMA_T$n[RMA_T$K==1.024],lag.max = 50)
# No periodicity

## In log-scale

par(mfrow=c(3,3))
spectrum(log(RMA$n[RMA$K==0.128]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==0.128]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==0.128]), method = "pgram")

spectrum(log(RMA$n[RMA$K==0.520]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==0.520]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==0.520]), method = "pgram")

spectrum(log(RMA$n[RMA$K==1.024]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==1.024]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==1.024]), method = "pgram")

## 

par(mfrow=c(3,3))
spec<-spectrum(log(RMA$n[RMA$K==0.128]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.128]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.128]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)

spec<-spectrum(log(RMA$n[RMA$K==0.520]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.520]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.520]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA$n[RMA$K==1.024]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==1.024]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==1.024]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)

## Perhaps I should zoom in the lower frequencies, with a max of 0.2=1/5
## With a sensible smoother
par(mfrow=c(3,3))
spec<-spectrum(log(RMA$n[RMA$K==0.128]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.128]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.128]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA$n[RMA$K==0.520]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.520]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.520]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA$n[RMA$K==1.024]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_I$n[RMA_I$K==1.024]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
spec<-spectrum(log(RMA_T$n[RMA_T$K==1.024]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.2))
# 1/0.02 is 50 so we should see long-period cycles


## Different method -- but we don't have enough points it seems

par(mfrow=c(3,3))
spectrum(RMA$n[RMA$K==0.128], method = "ar")
spectrum(RMA_I$n[RMA_I$K==0.128], method = "ar")
spectrum(RMA_T$n[RMA_T$K==0.128], method = "ar")

spectrum(RMA$n[RMA$K==0.520], method = "ar")
spectrum(RMA_I$n[RMA_I$K==0.520], method = "ar")
spectrum(RMA_T$n[RMA_T$K==0.520], method = "ar")

spectrum(RMA$n[RMA$K==1.024], method = "ar")
spectrum(RMA_I$n[RMA_I$K==1.024], method = "ar")
spectrum(RMA_T$n[RMA_T$K==1.024], method = "ar")

## Log-scale might be better -- but not really

par(mfrow=c(3,3))
spectrum(log(RMA$n[RMA$K==0.128]), method = "ar")
spectrum(log(RMA_I$n[RMA_I$K==0.128]), method = "ar")
spectrum(log(RMA_T$n[RMA_T$K==0.128]), method = "ar")

spectrum(log(RMA$n[RMA$K==0.520]), method = "ar")
spectrum(log(RMA_I$n[RMA_I$K==0.520]), method = "ar")
spectrum(log(RMA_T$n[RMA_T$K==0.520]), method = "ar")

spectrum(log(RMA$n[RMA$K==1.024]), method = "ar")
spectrum(log(RMA_I$n[RMA_I$K==1.024]), method = "ar")
spectrum(log(RMA_T$n[RMA_T$K==1.024]), method = "ar")