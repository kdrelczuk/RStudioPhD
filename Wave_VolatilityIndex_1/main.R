library(Quandl)
eurusd = Quandl("QUANDL/EURUSD", start_date = "2009-01-01", end_date = "2014-01-01")
eurusdRet = Quandl("QUANDL/EURUSD", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")
sandp = Quandl("YAHOO/INDEX_GSPC", start_date = "2009-01-01", end_date = "2014-01-01")
sandpRet = Quandl("YAHOO/INDEX_GSPC", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")
nasdaq = Quandl("YAHOO/INDEX_NDX", start_date = "2009-01-01", end_date = "2014-01-01")
nasdaqRet = Quandl("YAHOO/INDEX_NDX", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")

library(zoo)
rollapply(c(1,2,3,4,5,6,7,8,9,10,11),3,sum,by=3)

library(TFX)
QueryTrueFX(digits.secs=3)

library(wavethresh)

y = c(1,1,7,9,2,8,8,6)

ywd = wsd(eurusd$Rate[0:512])
ywd

names(ywd)

ywd$filter

library(wavelets)

head(eurusd$open[0:255])

v <- DJ.EX()$heavi
vnoise <- v + rnorm(length(v), 0, sd=sqrt(var(v))/6)
vnwst <- wst(vnoise)
plot(vnwst)
nrow(vnoise)
vnoise

require(biwavelet)
# Generate a synthetic time series 's' with the same power at three distinct periods
t1=sin(seq(from=0, to=2*5*pi, length=1000))
t2=sin(seq(from=0, to=2*15*pi, length=1000))
t3=sin(seq(from=0, to=2*40*pi, length=1000))
s=t1+t2+t3
# Compare non-corrected vs. corrected wavelet spectrum
wt1=wt(cbind(1:1000, s))
par(mfrow=c(1,2))
plot(wt1, type="power.corr.norm", main="Bias-corrected")
plot(wt1, type="power.norm", main="Not-corrected")
x1=xwt(cbind(1:1000, s), cbind(1:1000, s))
par(mfrow=c(1,2))
plot(x1, type="power.corr.norm", main="Bias-corrected")
plot(x1, type="power.norm", main="Not-corrected")

# Perform wavelet clustering of three time series
t1=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)))
t2=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)+0.1*pi))
t3=cbind(1:100, rnorm(100))
# Compute wavelet spectra
wt.t1=wt(t1)
wt.t2=wt(t2)
wt.t3=wt(t3)
# Store all wavelet spectra into array
w.arr=array(NA, dim=c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
w.arr[1, , ]=wt.t1$wave
w.arr[2, , ]=wt.t2$wave
w.arr[3, , ]=wt.t3$wave

w.arr.dis=wclust(w.arr)
plot(hclust(w.arr.dis$dist.mat, method="ward"), sub="", main="",
     ylab="Dissimilarity", hang=-1)



library(wmtsa)
## calculate the DWT of linear chirp 
linchirp <- make.signal("linchirp", n=1024)
eurusd$Rate[0:512]

y = c(1,1,7,9,2,8,8,6)

result   <- wavDWT(y, wavelet="haar", keep.series=TRUE)

names(result)
result$data
## plot the transform shifted for approximate zero 
## phase alignment 
plot(wavShift(result))

## plot summary 
eda.plot(result)

## summarize the transform 
summary(result)

wavStackPlot.default(wavMODWT(sunspots)$data)


eurDWT = wavDWT(eurusd$Rate[0:512],wavelet="d4", keep.series=TRUE)
eurDWT
plot(eurDWT)
eurMRD = wavMRD(eurDWT)
plot(eurMRD)
boxplot(eurMRD)

eurRetDWT = wavDWT(eurusd$Rate[512:1024],wavelet="d4")
eurRetMRD = wavMRD(eurRetDWT)
plot(eurRetMRD)
boxplot(as.matrix(eurRetMRD)[,1:9])

eurRetDWT = wavDWT(eurusdRet$Rate[0:512],wavelet="d4")
eurRetMRD = wavMRD(eurRetDWT)
plot(eurRetMRD)
boxplot(as.matrix(eurRetMRD)[,1:9])

plot(eurRetDWT)

boxplot(eurRetMRD)
boxplot(eurRetMRD[1:9])
boxplot(as.matrix(eurRetMRD)[,1:9])

test = c(10,10,10,10,10,10,198,10,10,10,10,10,10,10)
testDWT = wavDWT(test,wavelet="d2", keep.series=TRUE)
testMRD = wavMRD(testDWT)
par(mfrow=c(1,1))
plot(testMRD)
plot(testDWT)
boxplot(testMRD)

##############################################################################

library(Hmisc)
library(plyr)

sandp = Quandl("YAHOO/INDEX_GSPC", start_date = "2009-01-01", end_date = "2014-01-01")
nasdaq = Quandl("YAHOO/INDEX_NDX", start_date = "2009-01-01", end_date = "2014-01-01")
eurusd = Quandl("QUANDL/EURUSD", start_date = "2009-01-01", end_date = "2014-01-01")

sandpRet = Quandl("YAHOO/INDEX_GSPC", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")
nasdaqRet = Quandl("YAHOO/INDEX_NDX", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")
eurusdRet = Quandl("QUANDL/EURUSD", start_date = "2009-01-01", end_date = "2014-01-01",transformation="rdiff")

ticks = rev(sandp$Open)
#ticks = rev(nasdaq$Open)
#ticks = rev(eurusd$Rate)

#ticks = rev(sandpRet$Open)
#ticks = rev(nasdaqRet$Open)
#ticks = rev(eurusdRet$Rate)

getValIndex = function(data)
{
  dwt = wavDWT(data,wavelet="d4")
  mrd = wavMRD(dwt)
  #s = summary(mrd)
  sd(apply(mrd[,1:5],1,sum))
  #s$crystal.stat[1:5,7]
}

hasAnomaly = function(data)
{
  dwt = wavDWT(data,wavelet="d4")
  mrd = wavMRD(dwt)
  mrdMat = as.matrix(mrd)
  sdval1 = sd(mrdMat[,1])
  sdval2 = sd(mrdMat[,2])
  sdval3 = sd(mrdMat[,3])
  max(c(sum(mrdMat[,1]>sdval1*2),sum(mrdMat[,2]>sdval2*2),sum(mrdMat[,3]>sdval3*2)))
}

getRet = function(data)
{
  
  (data[length(data)]-data[1])/data[1]   
}

risk = rollapply(ticks,32,getRet,by=32)
ft = rollapply(ticks,32,getValIndex,by=32)
#colnames(ft) = c("D1","D2","D3","D4","D5")

all = cbind(ft,risk,abs(risk-Lag(risk,-1)))
colnames(all) = c("D","Lag0", "CH")
pairs(all)
cor(all,use="pairwise.complete.obs")


##############################################################################

mrdMat
apply(mrdMat[,1:5],1,sum)

hasAnomaly = function(data)
{
  dwt = wavDWT(data,wavelet="d4")
  mrd = wavMRD(dwt)
  mrdMat = as.matrix(mrd)
  sdval1 = sd(mrdMat[,1])
  sdval2 = sd(mrdMat[,2])
  sdval3 = sd(mrdMat[,3])
  c(sum(mrdMat[,1]>sdval1*2),sum(mrdMat[,2]>sdval2*2),sum(mrdMat[,3]>sdval3*2))
}

ft = rollapply(ticks,32,hasAnomaly,by=32)



x = c(1,2,3,4,5,6)
x[1:length(x)]
boxplot(as.matrix(sandpMRD)[,1:5])


getValIndex(sandp$Open[0:32])
s$crystal.stat[,7]

plot(sandpMRD)
boxplot(as.matrix(sandpMRD)[,1:5])

sandpDWT = wavDWT(sandp$Open[512:(512+32)],wavelet="d4")
sandpMRD = wavMRD(sandpDWT)
plot(sandpMRD)
boxplot(as.matrix(sandpMRD)[,1:5])

sandpDWT2 = wavDWT(sandp$Open[512:1024],wavelet="d4")
sandpMRD2 = wavMRD(sandpDWT2)
plot(sandpMRD2)
boxplot(as.matrix(sandpMRD2)[,1:9])

eda.plot(sandpDWT)
summary(sandpMRD)