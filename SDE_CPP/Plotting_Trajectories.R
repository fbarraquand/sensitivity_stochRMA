### R file for plotting trajectories of SDEs -- FB 18/04/2020

RMA = read.table("RMA.txt")
ncol(RMA)
names(RMA) = c("n","p","n2","p2","time","K")

RMA_I = read.table("RMA_I.txt")
ncol(RMA_I)
names(RMA_I) = c("n","p","n2","p2","time","K")

RMA_T = read.table("RMA_T.txt")
ncol(RMA_T)
names(RMA_T) = c("n","p","n2","p2","time","K")

pdf(file = "RMA_BifDiagK.pdf",width = 10, height = 10)
#png(file = "RMA_BifDiagK.png",width = 20, height = 20, units = "cm", res = 300)
par(mfrow=c(3,2),pch = 20, mar = c(4.2, 4 , 2, 1))

logn = log(RMA$n)
plot(RMA$K,log(RMA$n),col="black", type="p",xlim = c(min(RMA$K),max(RMA$K)),ylim = c(-45,max(logn)),
     xlab ="",ylab = "log(prey density)", main = "Holling") #xaxt = "n" to remove the axis ticks
lines(RMA$K,log(RMA$n2),col="red", type="p")

plot(RMA$K,RMA$n,col="black", type="p",
     xlab ="",ylab = "prey density", main = "Holling")
lines(RMA$K,RMA$n2,col="red", type="p")

plot(RMA_I$K,log(RMA_I$n),col="black", type="p",xlim = c(min(RMA$K),max(RMA$K)),ylim = c(-45,max(logn)),
    xlab ="",ylab = "log(prey density)", main = "Ivlev")
lines(RMA_I$K,log(RMA_I$n2),col="red", type="p")

plot(RMA_I$K,RMA_I$n,col="black", type="p",
     xlab ="",ylab = "prey density", main = "Ivlev")
lines(RMA_I$K,RMA_I$n2,col="red", type="p")

plot(RMA_T$K,log(RMA_T$n),col="black", type="p",xlim = c(min(RMA$K),max(RMA$K)),ylim = c(-45,max(logn)),
     xlab = "K", ylab = "log(prey density)", main = "tanh")
lines(RMA_T$K,log(RMA_T$n2),col="red", type="p")

plot(RMA_T$K,RMA_T$n,col="black", type="p",
     xlab = "K", ylab = "prey density", main = "tanh")
lines(RMA_T$K,RMA_T$n2,col="red", type="p")

dev.off()