### R file for plotting the power spectrum - FBarraquand 
### Inspiration for smoothing and labelling https://rstudio-pubs-static.s3.amazonaws.com/9428_1197bd003ebd43c49b429f22ea4f36e5.html

RMA = read.csv("timeseries_H.csv")
ncol(RMA)
names(RMA) = c("time","n","p","K")

RMA_I = read.csv("timeseries_I.csv")
ncol(RMA_I)
names(RMA_I) = c("time","n","p","K")

RMA_T = read.csv("timeseries_T.csv")
ncol(RMA_T)
names(RMA_T) = c("time","n","p","K")

plot(RMA$time[RMA$K==0.1],RMA$n[RMA$K==0.1],type="o",xlab = "Time",ylab = "Prey abundance")

### Zoom on the time series
par(mfrow=c(3,3))
plot(RMA$time[RMA$K==0.1&(RMA$time>100&RMA$time<250)],RMA$n[RMA$K==0.1&(RMA$time>100&RMA$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_I$time[RMA_I$K==0.1&(RMA_I$time>100&RMA_I$time<250)],RMA_I$n[RMA_I$K==0.1&(RMA_I$time>100&RMA_I$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_T$time[RMA_T$K==0.1&(RMA_T$time>100&RMA_T$time<250)],RMA_T$n[RMA_T$K==0.1&(RMA_T$time>100&RMA_T$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")

plot(RMA$time[RMA$K==0.5&(RMA$time>100&RMA$time<250)],RMA$n[RMA$K==0.5&(RMA$time>100&RMA$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_I$time[RMA_I$K==0.5&(RMA_I$time>100&RMA_I$time<250)],RMA_I$n[RMA_I$K==0.5&(RMA_I$time>100&RMA_I$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_T$time[RMA_T$K==0.5&(RMA_T$time>100&RMA_T$time<250)],RMA_T$n[RMA_T$K==0.5&(RMA_T$time>100&RMA_T$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")

plot(RMA$time[RMA$K==1&(RMA$time>100&RMA$time<250)],RMA$n[RMA$K==1&(RMA$time>100&RMA$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_I$time[RMA_I$K==1&(RMA_I$time>100&RMA_I$time<250)],RMA_I$n[RMA_I$K==1&(RMA_I$time>100&RMA_I$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")
plot(RMA_T$time[RMA_T$K==1&(RMA_T$time>100&RMA_T$time<250)],RMA_T$n[RMA_T$K==1&(RMA_T$time>100&RMA_T$time<250)],type="o",xlab = "Time",ylab = "Prey abundance")



#pdf(file = "spectrum.pdf",width=8,height=8)
par(mfrow=c(3,3))
spectrum(RMA$n[RMA$K==0.1],main = "Holling (K = 0.1)")
spectrum(RMA_I$n[RMA_I$K==0.1],main = "Ivlev (K = 0.1)")
spectrum(RMA_T$n[RMA_T$K==0.1],main = "tanh (K = 0.1)")

spectrum(RMA$n[RMA$K==0.5], main = "Holling (K = 0.5)")
spectrum(RMA_I$n[RMA_I$K==0.5], main = "Ivlev (K = 0.5)")
spectrum(RMA_T$n[RMA_T$K==0.5], main = "tanh (K = 0.5)")

spectrum(RMA$n[RMA$K==1], main = "Holling (K = 1)")
spectrum(RMA_I$n[RMA_I$K==1], main = "Ivlev (K = 1)")
spectrum(RMA_T$n[RMA_T$K==1], main = "tanh (K = 1)")
#dev.off()

#pdf(file = "spectrum_dB.pdf",width=8,height=8)
par(mfrow=c(3,3))
spectrum(RMA$n[RMA$K==0.1],log="dB",main = "Holling (K = 0.1)")
spectrum(RMA_I$n[RMA_I$K==0.1],log="dB",main = "Ivlev (K = 0.1)")
spectrum(RMA_T$n[RMA_T$K==0.1],log="dB",main = "tanh (K = 0.1)")

spectrum(RMA$n[RMA$K==0.5],log="dB", main = "Holling (K = 0.5)")
spectrum(RMA_I$n[RMA_I$K==0.5],log="dB", main = "Ivlev (K = 0.5)")
spectrum(RMA_T$n[RMA_T$K==0.5],log="dB", main = "tanh (K = 0.5)")

spectrum(RMA$n[RMA$K==1],log="dB", main = "Holling (K = 1)")
spectrum(RMA_I$n[RMA_I$K==1],log="dB", main = "Ivlev (K = 1)")
spectrum(RMA_T$n[RMA_T$K==1],log="dB", main = "tanh (K = 1)")
#dev.off()

#pdf(file = "ACF.pdf",width=8,height=8)
par(mfrow=c(3,3))
acf(RMA$n[RMA$K==0.1],lag.max=10000,main = "Holling (K = 0.1)")
acf(RMA_I$n[RMA_I$K==0.1],lag.max=10000,main = "Ivlev (K = 0.1)")
acf(RMA_T$n[RMA_T$K==0.1],lag.max=10000,main = "tanh (K = 0.1)")

acf(RMA$n[RMA$K==0.5],lag.max=10000, main = "Holling (K = 0.5)")
acf(RMA_I$n[RMA_I$K==0.5],lag.max=10000, main = "Ivlev (K = 0.5)")
acf(RMA_T$n[RMA_T$K==0.5],lag.max=10000, main = "tanh (K = 0.5)")

acf(RMA$n[RMA$K==1],lag.max=10000, main = "Holling (K = 1)", xlab =" Time lag (true units x 20)")
acf(RMA_I$n[RMA_I$K==1],lag.max=10000, main = "Ivlev (K = 1)", xlab =" Time lag (true units x 20)")
acf(RMA_T$n[RMA_T$K==1],lag.max=10000, main = "tanh (K = 1)", xlab =" Time lag (true units x 20)")
#dev.off()

#pdf(file = "ACF.pdf",width=8,height=8)
par(mfrow=c(3,3))
acf(RMA$n[RMA$K==0.1],lag.max=1000,main = "Holling (K = 0.1)")
acf(RMA_I$n[RMA_I$K==0.1],lag.max=1000,main = "Ivlev (K = 0.1)")
acf(RMA_T$n[RMA_T$K==0.1],lag.max=1000,main = "tanh (K = 0.1)")

acf(RMA$n[RMA$K==0.5],lag.max=1000, main = "Holling (K = 0.5)")
acf(RMA_I$n[RMA_I$K==0.5],lag.max=1000, main = "Ivlev (K = 0.5)")
acf(RMA_T$n[RMA_T$K==0.5],lag.max=1000, main = "tanh (K = 0.5)")

acf(RMA$n[RMA$K==1],lag.max=1000, main = "Holling (K = 1)", xlab =" Time lag (true units x 20)")
acf(RMA_I$n[RMA_I$K==1],lag.max=1000, main = "Ivlev (K = 1)", xlab =" Time lag (true units x 20)")
acf(RMA_T$n[RMA_T$K==1],lag.max=1000, main = "tanh (K = 1)", xlab =" Time lag (true units x 20)")
#dev.off()

## In log-scale (the spectrum itself, not the representation of the spectrum)

par(mfrow=c(3,3))
spectrum(log(RMA$n[RMA$K==0.1]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==0.1]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==0.1]), method = "pgram")

spectrum(log(RMA$n[RMA$K==0.5]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==0.5]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==0.5]), method = "pgram")

spectrum(log(RMA$n[RMA$K==1]), method = "pgram")
spectrum(log(RMA_I$n[RMA_I$K==1]), method = "pgram")
spectrum(log(RMA_T$n[RMA_T$K==1]), method = "pgram")
#Hardly changes anything


## Zoom of the periodogram on the interesting frequencies
pdf(file = "spectrum.pdf",width=8,height=8)
par(mfrow=c(3,3))
# https://stat.ethz.ch/pipermail/r-help/2006-January/086645.html
spectrum_n = spectrum(RMA$n[RMA$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 0.1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 0.1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 0.1)",ylab="Spectral density",xlab= "Frequency")

spectrum_n = spectrum(RMA$n[RMA$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 0.5)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 0.5)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 0.5)",ylab="Spectral density",xlab= "Frequency")

spectrum_n = spectrum(RMA$n[RMA$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 1)",ylab="Spectral density",xlab= "Frequency")
dev.off()

## Back-calculate the period used in 
1/min(spectrum_n$freq[spectrum_n$freq<0.01]) # max period in index units, 0.05 each,
1/max(spectrum_n$freq[spectrum_n$freq<0.01]) # min period in index units, 0.05 each

0.05/min(spectrum_n$freq[spectrum_n$freq<0.01]) # max period in original time units
0.05/max(spectrum_n$freq[spectrum_n$freq<0.01]) # min period in original time units
# frequency is below a certain level so periodicity is above a certain level (5 time units here)
# 0.01 in index unit frequency is 100 in index unit period = 100*0.05 = 5 in original time units
# 0.002 the dominant frequency is 500 in index unit period = 500*0.05 = 25 time units

## Logged representation

## Zoom of the periodogram on the interesting frequencies
par(mfrow=c(3,3))
spectrum_n = spectrum(RMA$n[RMA$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Holling (K = 0.1)")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Ivlev (K = 0.1)")
spectrum_n =spectrum(RMA_T$n[RMA_T$K==0.1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "tanh (K = 0.1)")

spectrum_n = spectrum(RMA$n[RMA$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Holling (K = 0.5)")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Ivlev (K = 0.5)")
spectrum_n =spectrum(RMA_T$n[RMA_T$K==0.5], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "tanh (K = 0.5)")

spectrum_n = spectrum(RMA$n[RMA$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Holling (K = 1)")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "Ivlev (K = 1)")
spectrum_n =spectrum(RMA_T$n[RMA_T$K==1], plot = FALSE)
plot(spectrum_n$freq[spectrum_n$freq<0.01],log(spectrum_n$spec[spectrum_n$freq<0.01]),type="l",main = "tanh (K = 1)")
## Log-representation kinda obscures everything for once ? 

## Smoothed spectrum
pdf(file = "spectrum_smoothed_and_zoomed.pdf",width=8,height=8)
par(mfrow=c(3,3))
spectrum_n = spectrum(RMA$n[RMA$K==0.1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 0.1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 0.1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==0.1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 0.1)",ylab="Spectral density",xlab= "Frequency")

spectrum_n = spectrum(RMA$n[RMA$K==0.5], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 0.5)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==0.5], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 0.5)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==0.5], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 0.5)",ylab="Spectral density",xlab= "Frequency")

spectrum_n = spectrum(RMA$n[RMA$K==1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Holling (K = 1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_I$n[RMA_I$K==1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "Ivlev (K = 1)",ylab="Spectral density",xlab= "Frequency")
spectrum_n = spectrum(RMA_T$n[RMA_T$K==1], plot = FALSE, spans = c(11,11,11))
plot(spectrum_n$freq[spectrum_n$freq<0.01]/0.05,2*spectrum_n$spec[spectrum_n$freq<0.01],type="l",main = "tanh (K = 1)",ylab="Spectral density",xlab= "Frequency")
dev.off()

## Usual smoothing -- not ideal

par(mfrow=c(3,3))
spec<-spectrum(log(RMA$n[RMA$K==0.1]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.1]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.1]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)

spec<-spectrum(log(RMA$n[RMA$K==0.5]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.5]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.5]), method = "pgram", kernel("daniell", c(40,20)),
         plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA$n[RMA$K==1]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_I$n[RMA_I$K==1]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)
spec<-spectrum(log(RMA_T$n[RMA_T$K==1]), method = "pgram", kernel("daniell", c(40,20)),
               plot = FALSE)
plot(spec)

## Perhaps I should zoom in the lower frequencies, with a max of 0.2=1/5
## With a sensible smoother
par(mfrow=c(3,3))
spec<-spectrum(log(RMA$n[RMA$K==0.1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Holling (K = 0.1)")
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Ivlev (K = 0.1)")
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "tanh (K = 0.1)")
spec<-spectrum(log(RMA$n[RMA$K==0.5]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Holling (K = 0.5)")
spec<-spectrum(log(RMA_I$n[RMA_I$K==0.5]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Ivlev (K = 0.5)")
spec<-spectrum(log(RMA_T$n[RMA_T$K==0.5]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "tanh (K = 0.5)")
spec<-spectrum(log(RMA$n[RMA$K==1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Holling (K = 1)")
spec<-spectrum(log(RMA_I$n[RMA_I$K==1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "Ivlev (K = 1)")
spec<-spectrum(log(RMA_T$n[RMA_T$K==1]), method = "pgram", kernel("daniell", c(50,30)),
              plot = FALSE)
plot(spec,xlim=c(0,0.01),main = "tanh (K = 1)")
# 1/0.02 is 50 so we should see long-period cycles


## Different method -- now we probably have enough points
par(mfrow=c(3,3))
spectrum(RMA$n[RMA$K==0.1], method = "ar")
spectrum(RMA_I$n[RMA_I$K==0.1], method = "ar")
spectrum(RMA_T$n[RMA_T$K==0.1], method = "ar")

spectrum(RMA$n[RMA$K==0.5], method = "ar")
spectrum(RMA_I$n[RMA_I$K==0.5], method = "ar")
spectrum(RMA_T$n[RMA_T$K==0.5], method = "ar")

spectrum(RMA$n[RMA$K==1], method = "ar")
spectrum(RMA_I$n[RMA_I$K==1], method = "ar")
spectrum(RMA_T$n[RMA_T$K==1], method = "ar")
# suggests it is predominantly red



