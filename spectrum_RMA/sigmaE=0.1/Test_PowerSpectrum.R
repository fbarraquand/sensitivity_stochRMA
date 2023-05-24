### Testing the spectrum reconstruction with artificial signals, FB 20/05/2022

## Using same expected periocity and TS length
size_vec = 40000
time_max = 2000
increment = time_max/size_vec
period = 20
period_index = round(period/increment)

x=seq(0.01,time_max,by = increment)
y = cos(2*pi*x/period) + 1*rnorm(size_vec)
plot(y[1:1000])

acf(y,lag.max=1000) #period is 20, so 400 in index units
spectrum(y) #zoom needed although there's a spike
spectrum_y = spectrum(y,spans=c(3,3)) # In that case spectrum(y) also works 
## We want period between 0 and 1000 or 2000 -- so frequency >1/(1000)
spectrum_y$spec[spectrum_y$freq>(1/1000)] # ? 
## or should we just zoom in the lower frequencies? 1/0.01 = 100 in index units = 5 in original units
plot(spectrum_y$spec[spectrum_y$freq<0.01])
plot(spectrum_y$spec[spectrum_y$freq<0.01],type="l")
plot(log(spectrum_y$spec[spectrum_y$freq<0.01]),type="l") #kinda obfuscate things ? 

## let's try Welch
library(gsignal)
fs = 20
nfft = fs/0.005 # 1/0.005 = 200 is then the max period we consider. 
spectrum_y = pwelch(y, detrend="none",fs=fs,nfft=nfft,overlap = 0.1)
# 0.2 = 1/5 is the max frequency we're interested in since we want to see oscillations in ~ 20 time units or more
plot(spectrum_y$freq[spectrum_y$freq<0.2],spectrum_y$spec[spectrum_y$freq<0.2],type="l",ylab="Spectral density",xlab= "Frequency")
## doesn't quite seem to work as expected... 

spectrum_y = pwelch(y, detrend="none",fs=fs,window = 5000,nfft=nfft,overlap = 0.1)
plot(spectrum_y$freq[spectrum_y$freq<0.2],spectrum_y$spec[spectrum_y$freq<0.2],type="l",ylab="Spectral density",xlab= "Frequency")
## much better


## check with another library
# https://www.rdocumentation.org/packages/bspec/versions/1.6/topics/welchPSD
library(bspec)
spec1 <- welchPSD(sunspots, seglength=200) 
plot(spec1$frequency, spec1$power,type="l", xlim=c(0,1))

## let's try that on our data
spec1 <- welchPSD(as.ts(y), seglength=5000) 
plot(spec1$frequency, spec1$power,type="l", xlim=c(0,0.01))
## OK that's better -- dominant frequency in index units is 1/400

