## ----setup, include=FALSE-----------------------------------------------------
remove(list=ls())
knitr::knit_hooks$set(purl = knitr::hook_purl)
library("readxl")
library("fUnitRoots")
library("astsa")
library("tidyverse")
library("fImport")
library("fBasics")
library("rugarch")
library("highfrequency")
library("xts")
library("zoo")
library("stats")
library("chron")
library("forecast")
library("quantmod")
library("astsa")
library("readxl")
library("tidyverse")
library("timeSeries")
library("fImport")
library("fBasics")
library("rugarch")
library("aTSA")
library("stats")
library(sandwich)
library(dplyr)
library(tidyr)
library(PortfolioAnalytics)
library(DEoptim)
library(ROI)
require(ROI.plugin.quadprog)
require(ROI.plugin.glpk)

## -----------------------------------------------------------------------------
getSymbols('VT', from='2008-06-26', to='2023-03-31',warnings=FALSE, 
           auto.assign=TRUE,  periodicity='daily')
VTd = `VT` $ `VT.Adjusted`
VTd = na.omit(VTd)

getSymbols('EEM', from='2003-04-14', to='2023-03-31',warnings=FALSE, 
           auto.assign=TRUE,  periodicity='daily')
EEMd = `EEM` $ `EEM.Adjusted`
EEMd = na.omit(EEMd)

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(2,3))
plot(VTd, main='VT - DAILY PRICES')
plot(log(VTd), main='VT - DAILY LOG PRICES')
plot(diff(log(VTd)), main='VT - DAILY LOG RETURN')
acf(VTd, lag=50)
acf(log(VTd), lag=50)
acf(na.omit(diff(log(VTd))), lag=50)

## ----warning=FALSE------------------------------------------------------------
dlvts = na.omit(diff(log(VTd)))
outB = Box.test(dlvts,1:20)
round(outB$p.value, 5)

## -----------------------------------------------------------------------------
par(mfrow=c(2,3))
dlvts = na.omit(diff(log(VTd)))
plot(density(dlvts), main='VTd- DAILY RETURNS DENSITY')
qqnorm(dlvts,main="Normal Q-Q Plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqline(dlvts,col="steelblue",lwd=2)
acf(dlvts,20,main='')
dlvts2=(dlvts)^2
plot(dlvts2, main='VTd - SQUARED LOG.RET')
acf(dlvts2, 20)
pacf(dlvts2, 20)      

## -----------------------------------------------------------------------------
out1 = adfTest(log(VTd), lag=21, type='c')
summary(out1@test$lm)

## -----------------------------------------------------------------------------
out1 = adfTest(log(VTd), lag=20, type='nc')
summary(out1@test$lm)

## ----include=FALSE------------------------------------------------------------
#important because the number of lags considered in an ADF (Augmented Dickey-Fuller test) represents the amount of lagging of the error term that is included in the model. In practice, each additional delay corresponds to a greater consideration of the past variations of the historical series in the statistical model used for the test.

## -----------------------------------------------------------------------------
adfTest(log(VTd), lags=20, type='nc')

## ----warning=FALSE------------------------------------------------------------
out1 = adfTest(dlvts, lag=19, type='nc')
summary(out1@test$lm)
adfTest(dlvts, lags=19, type='nc')

## -----------------------------------------------------------------------------
p_VT = log(window(VTd, end="2022-03-31")) #log-prices
p_VT = as.vector(p_VT$VT.Adjusted[,1])
T = length(p_VT)
r_VT=100*(p_VT[2:T]-p_VT[1:T-1]) #log returns &
par(mfrow=c(1,1))
tsplot(r_VT, ylab='VT') # plot log returns %

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
acf(r_VT,lag.max=40,xlab="",ylab="",main="VT ret")
acf(abs(r_VT),lag.max=40,xlab="",ylab="",main="VT abs")
acf(abs(r_VT)^2,lag.max=40,xlab="",ylab="",main="vt sqr")

## ----warning=FALSE------------------------------------------------------------
par(mfrow=c(1,1))
rollV=rollStats(as.timeSeries(r_VT),62,FUN=var)
tsplot(rollV,ylab="",xlab="")
rollV=rollStats(as.timeSeries(r_VT),126,FUN=var)
tsplot(rollV,ylab="",xlab="")

## -----------------------------------------------------------------------------
y2 = (r_VT - mean(r_VT))^2 #returns in deviation from mean,squared
Box.test(y2, lag = 5, type = "Ljung-Box")
Box.test(y2, lag = 10, type = "Ljung-Box")

## -----------------------------------------------------------------------------
# consider 1, 2, and 3 lags for the LM test
y2=as.timeSeries(y2)
y2L1 <- lag(y2,1)
y2L2 <- lag(y2,2)
y2L3 <- lag(y2,3)
# three test statistics and corresponding P-values
# 
# First lag:
out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
# using T-2 because some observations are lost in lags 1,2,3 and to compare I choose 2
pchisq(LMARCH1,1,lower.tail=FALSE) # test only on the upper tail
# Second lag:
out2 <- lm(y2 ~ y2L1 + y2L2)
sum2 <- summary(out2)
LMARCH2 <- sum2$r.squared*(T-2)
pchisq(LMARCH2,2,lower.tail=FALSE)
# Third lag:
out3 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum3 <- summary(out3)
LMARCH3 <- sum3$r.squared*(T-2)
pchisq(LMARCH3,3,lower.tail=FALSE)

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
tsplot(r_VT, ylab='VT') # plot of log percentage returns
spec0 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
# estimate model saving output
fit0 <- ugarchfit(spec0, r_VT)
# estimate model with output on screen
fit0@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(2,1))
tsplot(r_VT, type='l', ylab="", xlab="")
tsplot(fit0@fit$sigma, type='l',ylab="",xlab="")

## -----------------------------------------------------------------------------
sim0<-ugarchsim(fit0, n.sim = 5000)
par(mfrow=c(2,1))
tsplot(sim0@simulation$seriesSim,type="l",ylab="",xlab="")
tsplot(sim0@simulation$sigmaSim,type="l",ylab="",xlab="")

## -----------------------------------------------------------------------------
resst0 <- residuals(fit0,standardize=TRUE)
plot(resst0)
# analysis for GARCH effects
y2<-(as.timeSeries(resst0))^2
# consider 1, 2, and 3 lags for the LM test
y2L1 <- lag(y2,1)
y2L2 <- lag(y2,2)
y2L3 <- lag(y2,3)

# Evaluate the regression to understand if zt ^ 2 depends on its lags.
# If alpha = 0, zt is homoscedastic.
# If alpha is different from zero, zt is heteroscedastic.
# Test with 1,2,3 delays, the test statistic changes.
# The null hypothesis is that the lags are jointly zero.
# If the test is correctly specified, I expect to have
# no further evidence of ARCH effects.
# three test statistics and their P-values
out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,1,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,2,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,3,lower.tail=FALSE)

# In this case, we do not reject the null hypothesis, the coefficients are zero,
# in the case of the first, only one coefficient.

## -----------------------------------------------------------------------------
# Ljung-Box test (stats) at 5 and 10 lags
Box.test(as.numeric(y2), lag = 5, type = "Ljung-Box")
Box.test(as.numeric(y2), lag = 10, type = "Ljung-Box")

## -----------------------------------------------------------------------------
# Correlograms
par(mfrow=c(2,3))
acf(r_VT,lag.max=40,xlab="",ylab="",main="VT ret")
acf(abs(r_VT),lag.max=40,xlab="",ylab="",main="VT abs")
acf(abs(r_VT)^2,lag.max=40,xlab="",ylab="",main="VT sqr")
acf(as.numeric(resst0),lag.max=40,xlab="",ylab="",main="VT ret")
acf(abs(as.numeric(resst0)),lag.max=40,xlab="",ylab="",main="VT abs")
acf(abs(as.numeric(resst0))^2,lag.max=40,xlab="",ylab="",main="VT sqr")

## ----fig.align='center'-------------------------------------------------------
# Evaluating the density of residuals
# QQ-plot against alternative distributions
d0<-density(resst0)
# Residual density
par(mfrow=c(1,1))
plot (d0,main ="VT stdres density",ylab ="")
# QQ-plot against normal
qqnorm(resst0,main="Normal Q-Q Plot",xlab="Theoretical Quantiles",ylab ="Sample Quantiles")
qqline(resst0,col="steelblue",lwd=2)

## -----------------------------------------------------------------------------
y=r_VT-mean(r_VT)
nnr=y*(y<0)# per negative sign bias
ppr=y*(y>0)# per positive signi bias
dn=y<0# per sign bias
nnr=as.timeSeries(nnr)
ppr=as.timeSeries(ppr)
dn=as.timeSeries(dn)
nnr=lag(nnr,1)
ppr=lag(ppr,1)
dn=lag(dn,1)

out.sbt=lm(resst0^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resst0^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resst0^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resst0^2 ~ dn + nnr + ppr)
summary(out.jsbt)

## -----------------------------------------------------------------------------
p_VT = log(window(VTd, end="2022-03-31")) #log-prices
p_VT = as.vector(p_VT$VT.Adjusted[,1])
T = length(p_VT)
r_VT=100*(p_VT[2:T]-p_VT[1:T-1]) #log returns %
y = r_VT - mean(r_VT)

spec0n<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
mean.model = list(armaOrder = c(0, 0),
include.mean = TRUE), distribution.model="norm")
fit0n<-ugarchfit(spec0n,r_VT)
resst0n <- residuals(fit0n,standardize=TRUE)
fit0n@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
specgn <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
fitgn <- ugarchfit(specgn,r_VT)
resstgn <- residuals(fitgn,standardize=TRUE)
fitgn@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
specan <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),distribution.model="norm")
fitan <- ugarchfit(specan,r_VT)
resstan <- residuals(fitan,standardize=TRUE)
fitan@fit$robust.matcoef

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(resst0n)
plot(resstgn)
plot(resstan)
par(mfrow=c(1,1))

ks.test(as.matrix(resst0n),"pnorm",0,1)

## -----------------------------------------------------------------------------
spec0ta <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),distribution.model="sstd")
fit0ta <- ugarchfit(spec0ta,r_VT)
resst0ta <- residuals(fit0ta,standardize=TRUE)

specgta <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="sstd")
fitgta <- ugarchfit(specgta,r_VT)
resstgta <- residuals(fitgta,standardize=TRUE)

specata <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="sstd")
fitata <- ugarchfit(specata,r_VT)
resstata <- residuals(fitata,standardize=TRUE)

## -----------------------------------------------------------------------------
fit0ta@fit$robust.matcoef
fitgta@fit$robust.matcoef
fitata@fit$robust.matcoef

## -----------------------------------------------------------------------------
specata <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model="sstd")
fitata <- ugarchfit(specata,r_VT)
resstata <- residuals(fitata,standardize=TRUE)
fitata@fit$robust.matcoef

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(resst0ta)
plot(resstgta)
plot(resstata)
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
p=0.01*(1:99)
qemp=quantile(resst0ta,probs=p)
qteor=qdist(distribution="sstd",p,mu=0,sigma=1,skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6])
qqplot(qteor,qemp,xlab="Theoretical Quantiles",ylab ="Sample Quantiles")
qqline(qemp, distribution = function(p) 
   qdist(distribution="sstd",p,mu=0,sigma=1,
                  skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6]))

## -----------------------------------------------------------------------------
# H0: Skew-T
y=as.matrix(resst0ta)
y=as.matrix(sort(y))
ks.test(as.matrix(resst0ta),pdist(distribution="sstd",y,mu=0,sigma=1,lambda=-0.5,
skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6]))

## -----------------------------------------------------------------------------
# grafici varianze
par(mfrow=c(1,3))
sigma0n <- sigma(fit0n)
sigmagn <- sigma(fitgn)
sigmaan <- sigma(fitan)
par(mfrow=c(1,3))
plot(as.numeric(sigma0n),type="l",ylab="",xlab="")
title("GARCH")
plot(as.numeric(sigmagn),type="l",ylab="",xlab="")
title("GJR-GARCH(1,1)")
plot(as.numeric(sigmaan),type="l",ylab="",xlab="")
title("APARCH(1,1)")
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
y=r_VT-mean(r_VT)
nnr=y*(y<0)# for negative sign bias
ppr=y*(y>0)# for positive signi bias
dn=y<0# for sign bias
nnr=as.timeSeries(nnr)
ppr=as.timeSeries(ppr)
dn=as.timeSeries(dn)
nnr=lag(nnr,1)
ppr=lag(ppr,1)
dn=lag(dn,1)

out.sbt=lm(resst0ta^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resst0ta^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resst0ta^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resst0ta^2 ~ dn + nnr + ppr)
summary(out.jsbt)

out.sbt=lm(resstgta^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resstgta^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resstgta^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resstgta^2 ~ dn + nnr + ppr)
summary(out.jsbt)

out.sbt=lm(resstata^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resstata^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resstata^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resstata^2 ~ dn + nnr + ppr)
summary(out.jsbt)

## ----fig.align='center'-------------------------------------------------------
ICall<-cbind(infocriteria(fit0n),infocriteria(fitgn),infocriteria(fitan), infocriteria(fit0ta),infocriteria(fitgta),infocriteria(fitata))
colnames(ICall)<-c("GARCH N","GJR N","APARCH N","GARCH TA","GJR TA","APARCH TA")
ICall_df <- data.frame(ICall)
print(ICall_df)

## -----------------------------------------------------------------------------
p_VT = log(window(VTd, end="2023-03-31")) #log-prices
p_VT = as.vector(p_VT$VT.Adjusted[,1])
T = length(p_VT)
r_VT=100*(p_VT[2:T]-p_VT[1:T-1]) #log returns %
y = r_VT - mean(r_VT)

# GARCH(1,1) Normal
for0n = ugarchroll(spec0n,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
for0nd = as.data.frame(for0n)
s0nf = for0nd$Sigma

# GARCH(1,1) skew-t
for0ta = ugarchroll(spec0ta,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
for0tad = as.data.frame(for0ta)
s0taf = for0tad$Sigma

# GJRGARCH(1,1) Normal
forgn = ugarchroll(specgn,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forgnd = as.data.frame(forgn)
sgnf = forgnd$Sigma

# GJRGARCH(1,1) skew-t
forgta = ugarchroll(specgta,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forgtad = as.data.frame(forgta)
sgtaf = forgtad$Sigma

# APARCH(1,1) Normal
foran = ugarchroll(specan,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forand = as.data.frame(foran)
sanf = forand$Sigma

# APARCH(1,1) skew-t
forata = ugarchroll(specata,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
foratad = as.data.frame(forata)
sataf = foratad$Sigma

## ----fig.align='center'-------------------------------------------------------
fall = cbind(s0nf,s0taf,sgnf,sgtaf,sanf,sataf)
par(mfrow=c(1,1))
matplot(fall, type = c("b"),pch=1,col = 1:6)
f0 = cbind(s0nf,s0taf)
fg = cbind(sgnf,sgtaf)
fa = cbind(sanf,sataf)
matplot(f0, type=('b'),pch=1,col=1:2)
matplot(fg, type=c('b'),pch=1,col=1:2)
matplot(fa, type=c('b'),pch=1,col=1:2)

## -----------------------------------------------------------------------------
rf = for0nd$Realized

lN = (s0nf^2 -rf^2)^2 # GARCH
lTA= (s0taf^2 -rf^2)^2 # GARCH skew-t
lGN = (sgnf^2 -rf^2)^2 # GJRGARCH
lGTA = (sgtaf^2 -rf^2)^2 # GJRGARCH skew-t
lAN = (sanf^2 -rf^2)^2 # GJRGARCH t
lATA = (sataf^2 -rf^2)^2 # APARCH

## -----------------------------------------------------------------------------
dNTA = lN-lTA 
dNGN = lN-lGN 
dNGTA = lN-lGTA 
dNAN = lN-lAN 
dNATA = lN-lATA 

dTAGN = lTA-lGN 
dTAGTA = lTA-lGTA 
dTAAN = lTA-lAN 
dTAATA = lTA-lATA 

dGNGTA = lGN-lGTA 
dGNAN = lGN-lAN 
dGNATA = lGN-lATA 

dGTAAN = lGTA-lAN 
dGTAATA = lGTA-lATA 

dANATA = lAN-lATA

## -----------------------------------------------------------------------------
m<-floor(0.75*((NROW(rf))^(1/3)))

x1<-as.vector(matrix(1,nrow=NROW(rf)))

VNTA = NeweyWest(lm(dNTA~x1-1),lag=m,prewhite=0)
VNGN = NeweyWest(lm(dNGN~x1-1),lag=m,prewhite=0)
VNGTA = NeweyWest(lm(dNGTA~x1-1),lag=m,prewhite=0)
VNAN = NeweyWest(lm(dNAN~x1-1),lag=m,prewhite=0)
VNATA = NeweyWest(lm(dNATA~x1-1),lag=m,prewhite=0)

VTAGN = NeweyWest(lm(dTAGN~x1-1),lag=m,prewhite=0)
VTAGTA = NeweyWest(lm(dTAGTA~x1-1),lag=m,prewhite=0)
VTAAN = NeweyWest(lm(dTAAN~x1-1),lag=m,prewhite=0)
VTAATA = NeweyWest(lm(dTAATA~x1-1),lag=m,prewhite=0)

VGNGTA = NeweyWest(lm(dGNGTA~x1-1),lag=m,prewhite=0)
VGNAN = NeweyWest(lm(dGNAN~x1-1),lag=m,prewhite=0)
VGNATA = NeweyWest(lm(dGNATA~x1-1),lag=m,prewhite=0)

VGTAAN = NeweyWest(lm(dGTAAN~x1-1),lag=m,prewhite=0)
VGTAATA = NeweyWest(lm(dGTAATA~x1-1),lag=m,prewhite=0)

VANATA = NeweyWest(lm(dANATA~x1-1),lag=m,prewhite=0)

DM<-matrix(0,nrow=6,ncol=6)
# Output -> robust variances regressors

# loss = model in row minus model in column
colnames(DM)<-c("GARCHnorm","GARCH TA","GJRnorm","GJR ta","APARCHnorm","APARCH TA")
rownames(DM)<-c("GARCHnorm","GARCH TA","GJRnorm","GJR ta","APARCHnorm","APARCH TA")

DM[1, 2] <- mean(dNTA) / sqrt(VNTA)
DM[1, 3] <- mean(dNGN) / sqrt(VNGN)
DM[1, 4] <- mean(dNGTA) / sqrt(VNGTA)
DM[1, 5] <- mean(dNAN) / sqrt(VNAN)
DM[1, 6] <- mean(dNATA) / sqrt(VNATA)
DM[2, 3] <- mean(dTAGN) / sqrt(VTAGN)
DM[2, 4] <- mean(dTAGTA) / sqrt(VTAGTA)
DM[2, 5] <- mean(dTAAN) / sqrt(VTAAN)
DM[2, 6] <- mean(dTAATA) / sqrt(VTAATA)
DM[3, 4] <- mean(dGNGTA) / sqrt(VGNGTA)
DM[3, 5] <- mean(dGNAN) / sqrt(VGNAN)
DM[3, 6] <- mean(dGNATA) / sqrt(VGNATA)
DM[4, 5] <- mean(dGTAAN) / sqrt(VGTAAN)
DM[4, 6] <- mean(dGTAATA) / sqrt(VGTAATA)
DM[5, 6] <- mean(dANATA) / sqrt(VANATA)
DM

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(2,3))
plot(EEMd, main='EEMd - DAILY PRICES')
plot(log(EEMd), main='EEMd - DAILY LOG PRICES')
plot(diff(log(EEMd)), main='EEMd - DAILY LOG RETURN')
acf(EEMd, lag=50)
acf(log(EEMd), lag=50)
acf(na.omit(diff(log(EEMd))), lag=50)

## ----warning=FALSE------------------------------------------------------------
dlEEMs = na.omit(diff(log(EEMd)))
outB = Box.test(dlEEMs,1:20)
round(outB$p.value, 5)

## -----------------------------------------------------------------------------
par(mfrow=c(2,3))
dlEEMs = na.omit(diff(log(EEMd)))
plot(density(dlEEMs), main='EEMd- DAILY RETURNS DENSITY')
qqnorm(dlEEMs,main="Normal Q-Q Plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqline(dlEEMs,col="steelblue",lwd=2)
acf(dlEEMs,20,main='')
dlEEMs2=(dlEEMs)^2
plot(dlEEMs2, main='EEMd - SQUARED LOG.RET')
acf(dlEEMs2, 20)
pacf(dlEEMs2, 20)      

## -----------------------------------------------------------------------------
out1 = adfTest(log(EEMd), lag=21, type='c')
summary(out1@test$lm)

## -----------------------------------------------------------------------------
out1 = adfTest(log(EEMd), lag=18, type='c')
summary(out1@test$lm)

## ----include=FALSE------------------------------------------------------------
#important because the number of lags considered in an ADF (Augmented Dickey-Fuller test) represents the amount of lagging of the error term that is included in the model. In practice, each additional delay corresponds to a greater consideration of the past variations of the historical series in the statistical model used for the test.

## -----------------------------------------------------------------------------
adfTest(log(EEMd), lag=18, type='c')

## ----warning=FALSE------------------------------------------------------------
out1 = adfTest(dlEEMs, lag=17, type='nc')
summary(out1@test$lm)
adfTest(dlEEMs, lag=17, type='nc')

## -----------------------------------------------------------------------------
p_EEM = log(window(EEMd, end="2022-03-31")) #log-prices
p_EEM = as.vector(p_EEM$EEM.Adjusted[,1])
T = length(p_EEM)
r_EEM=100*(p_EEM[2:T]-p_EEM[1:T-1]) #log percentage returns
par(mfrow=c(1,1))
tsplot(r_EEM, ylab='EEM') # plot of log percentage returns

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
acf(r_EEM,lag.max=40,xlab="",ylab="",main="EEM ret")
acf(abs(r_EEM),lag.max=40,xlab="",ylab="",main="EEM abs")
acf(abs(r_EEM)^2,lag.max=40,xlab="",ylab="",main="EEM sqr")

## ----warning=FALSE------------------------------------------------------------
par(mfrow=c(1,1))
rollV=rollStats(as.timeSeries(r_EEM),62,FUN=var)
tsplot(rollV,ylab="",xlab="")
rollV=rollStats(as.timeSeries(r_EEM),126,FUN=var)
tsplot(rollV,ylab="",xlab="")

## -----------------------------------------------------------------------------
y2 = (r_EEM - mean(r_EEM))^2 # returns squared, deviation from the mean
Box.test(y2, lag = 5, type = "Ljung-Box")
Box.test(y2, lag = 10, type = "Ljung-Box")

## -----------------------------------------------------------------------------
# consider 1, 2, and 3 lags for the LM test
y2=as.timeSeries(y2)
y2L1 <- lag(y2,1)
y2L2 <- lag(y2,2)
y2L3 <- lag(y2,3)
# three test statistics and corresponding P-values
# 
# First lag:
out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
# using T-2 because we lose some observations in lags 1,2,3 and for comparison, I choose 2
pchisq(LMARCH1,1,lower.tail=FALSE)# test only on the upper tail
# Second lag:
out2 <- lm(y2 ~ y2L1 + y2L2)
sum2 <- summary(out2)
LMARCH2 <- sum2$r.squared*(T-2)
pchisq(LMARCH2,2,lower.tail=FALSE)
# Third lag:
out3 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum3 <- summary(out3)
LMARCH3 <- sum3$r.squared*(T-2)
pchisq(LMARCH3,3,lower.tail=FALSE)

## -----------------------------------------------------------------------------
spec0<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
# estimate model saving output
fit0<-ugarchfit(spec0,r_EEM)
# estimate model with output on screen
fit0@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(2,1))
tsplot(r_EEM, type='l', ylab="", xlab="")
tsplot(fit0@fit$sigma, type='l',ylab="",xlab="")

## -----------------------------------------------------------------------------
sim0<-ugarchsim(fit0, n.sim = 5000)
par(mfrow=c(2,1))
tsplot(sim0@simulation$seriesSim,type="l",ylab="",xlab="")
tsplot(sim0@simulation$sigmaSim,type="l",ylab="",xlab="")

## -----------------------------------------------------------------------------
resst0 <- residuals(fit0,standardize=TRUE)
plot(resst0)
# analysis for GARCH effects
y2 <- (as.timeSeries(resst0))^2
# consider 1, 2, and 3 lags for the LM test
y2L1 <- lag(y2,1)
y2L2 <- lag(y2,2)
y2L3 <- lag(y2,3)

# Evaluate the regression to understand if z_t^2 depends on its lags.
# If alpha = 0 then z_t is homoschedastic.
# If alpha is not 0 then z_t is heteroschedastic.
# Try with 1,2,3 lags, changing the test statistic.
# The null hypothesis is that the lags are jointly zero.
# If the test is correctly specified, I expect to have
# no further evidence of ARCH effects.
# three test statistics and corresponding P-values
out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,1,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,2,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,3,lower.tail=FALSE)

# In this case, we do not reject the null hypothesis, the coefficients are zero.

## -----------------------------------------------------------------------------
spec22<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
# estimate model saving output
fit22<-ugarchfit(spec22,r_EEM)
# estimate model with output on screen
fit22@fit$robust.matcoef

resst22 <- residuals(fit22,standardize=TRUE)
# analysis for GARCH effects
y2 <- (as.timeSeries(resst22))^2
# consider 1, 2, and 3 lags for the LM test
y2L1 <- lag(y2,1)
y2L2 <- lag(y2,2)
y2L3 <- lag(y2,3)

out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,1,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,2,lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1,3,lower.tail=FALSE)

## -----------------------------------------------------------------------------
spec22<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
# estimate model saving output
fit22<-ugarchfit(spec22, r_EEM)
# estimate model with output on screen
fit22@fit$robust.matcoef

resst22 <- residuals(fit22, standardize=TRUE)
# analysis for GARCH effects
y2 <- (as.timeSeries(resst22))^2
# consider 1, 2, and 3 lags for the LM test
y2L1 <- lag(y2, 1)
y2L2 <- lag(y2, 2)
y2L3 <- lag(y2, 3)

out1 <- lm(y2 ~ y2L1)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1, 1, lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1, 2, lower.tail=FALSE)
out1 <- lm(y2 ~ y2L1 + y2L2 + y2L3)
sum1 <- summary(out1)
LMARCH1 <- sum1$r.squared*(T-2)
pchisq(LMARCH1, 3, lower.tail=FALSE)

## -----------------------------------------------------------------------------
y2 <- (as.timeSeries(resst0))^2
# Ljung-Box test (stats) at 5 and 10 lags
Box.test(as.numeric(y2), lag = 5, type = "Ljung-Box")
Box.test(as.numeric(y2), lag = 10, type = "Ljung-Box")

## ----fig.align='center'-------------------------------------------------------
# correlograms
par(mfrow=c(2,3))
acf(r_EEM, lag.max=40, xlab="", ylab="", main="EEM ret")
acf(abs(r_EEM), lag.max=40, xlab="", ylab="", main="EEM abs")
acf(abs(r_EEM)^2, lag.max=40, xlab="", ylab="", main="EEM sqr")
acf(as.numeric(resst0), lag.max=40, xlab="", ylab="", main="EEM ret")
acf(abs(as.numeric(resst0)), lag.max=40, xlab="", ylab="", main="EEM abs")
acf(abs(as.numeric(resst0))^2, lag.max=40, xlab="", ylab="", main="EEM sqr")

## ----fig.align='center'-------------------------------------------------------
# Evaluation of residual density
# QQ-plot against alternative distributions
d0 <- density(resst0)
# Residual density
par(mfrow=c(1,1))
plot(d0, main="EEM stdres density", ylab="")
# Thick tails, density of standardized residuals
# QQ-plot against normal distribution
qqnorm(resst0, main="Normal Q-Q Plot", xlab="Theoretical Quantiles", ylab="Sample Quantiles")
qqline(resst0, col="steelblue", lwd=2)

## -----------------------------------------------------------------------------
y = r_EEM - mean(r_EEM)
nnr = y * (y < 0) # for negative sign bias
ppr = y * (y > 0) # for positive sign bias
dn = y < 0 # for sign bias
nnr = as.timeSeries(nnr)
ppr = as.timeSeries(ppr)
dn = as.timeSeries(dn)
nnr = lag(nnr, 1)
ppr = lag(ppr, 1)
dn = lag(dn, 1)

out.sbt = lm(resst0^2 ~ dn)
summary(out.sbt)
out.nsbt = lm(resst0^2 ~ nnr)
summary(out.nsbt)
out.psbt = lm(resst0^2 ~ ppr)
summary(out.psbt)
out.jsbt = lm(resst0^2 ~ dn + nnr + ppr)
summary(out.jsbt)

## -----------------------------------------------------------------------------
p_EEM = log(window(EEMd, end="2022-03-31")) # log-prices
p_EEM = as.vector(p_EEM$EEM.Adjusted[,1])
T = length(p_EEM)
r_EEM = 100 * (p_EEM[2:T] - p_EEM[1:T-1]) # log percentage returns
y = r_EEM - mean(r_EEM)

spec0n <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
fit0n <- ugarchfit(spec0n, r_EEM)
resst0n <- residuals(fit0n, standardize=TRUE)
fit0n@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
specgn <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
fitgn <- ugarchfit(specgn, r_EEM)
resstgn <- residuals(fitgn, standardize=TRUE)
fitgn@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
specgn <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model="norm")
fitgn <- ugarchfit(specgn, r_EEM)
resstgn <- residuals(fitgn, standardize=TRUE)
fitgn@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
specan <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="norm")
fitan <- ugarchfit(specan, r_EEM)
fitan@fit$robust.matcoef

## -----------------------------------------------------------------------------
specan <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model="norm")
fitan <- ugarchfit(specan, r_EEM)
resstan <- residuals(fitan, standardize=TRUE)
fitan@fit$robust.matcoef

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(resst0n)
plot(resstgn)
plot(resstan)
par(mfrow=c(1,1))

ks.test(as.matrix(resst0n),"pnorm",0,1)

## -----------------------------------------------------------------------------
spec0ta <- ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="sstd")
fit0ta <- ugarchfit(spec0ta,r_EEM)
resst0ta <- residuals(fit0ta,standardize=TRUE)

specgta <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="sstd")
fitgta <- ugarchfit(specgta,r_EEM)
resstgta <- residuals(fitgta,standardize=TRUE)

specata <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model="sstd")
fitata <- ugarchfit(specata,r_EEM)
resstata <- residuals(fitata,standardize=TRUE)

## -----------------------------------------------------------------------------
fit0ta@fit$robust.matcoef
fitgta@fit$robust.matcoef
fitata@fit$robust.matcoef

## -----------------------------------------------------------------------------
specgta <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model="sstd")
fitgta <- ugarchfit(specgta,r_EEM)
resstgta <- residuals(fitgta,standardize=TRUE)

specata <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model="sstd")
fitata <- ugarchfit(specata,r_EEM)
resstata <- residuals(fitata,standardize=TRUE)

## -----------------------------------------------------------------------------
fit0ta@fit$robust.matcoef
fitgta@fit$robust.matcoef
fitata@fit$robust.matcoef

## -----------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(resst0ta)
plot(resstgta)
plot(resstata)
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
p=0.01*(1:99)
qemp=quantile(resst0ta,probs=p)
qteor=qdist(distribution="sstd",p,mu=0,sigma=1,skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6])
qqplot(qteor,qemp,xlab="Theoretical Quantiles",ylab ="Sample Quantiles")
qqline(qemp, distribution = function(p) qdist(distribution="sstd",p,mu=0,sigma=1,skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6]))

## -----------------------------------------------------------------------------
# H0: Skew-T
y=as.matrix(resst0ta)
y=as.matrix(sort(y))
ks.test(as.matrix(resst0ta),pdist(distribution="sstd",y,mu=0,sigma=1,lambda=-0.5,
skew=fit0ta@fit$coef[5],shape=fit0ta@fit$coef[6]))

## -----------------------------------------------------------------------------
# variance plots
par(mfrow=c(1,3))
sigma0n <- sigma(fit0n)
sigmagn <- sigma(fitgn)
sigmaan <- sigma(fitan)
par(mfrow=c(1,3))
plot(as.numeric(sigma0n),type="l",ylab="",xlab="")
title("GARCH")
plot(as.numeric(sigmagn),type="l",ylab="",xlab="")
title("GJR-GARCH(1,1)")
plot(as.numeric(sigmaan),type="l",ylab="",xlab="")
title("APARCH(1,1)")
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
y=r_EEM-mean(r_EEM)
nnr=y*(y<0)# for negative sign bias
ppr=y*(y>0)# for positive signi bias
dn=y<0# for sign bias
nnr=as.timeSeries(nnr)
ppr=as.timeSeries(ppr)
dn=as.timeSeries(dn)
nnr=lag(nnr,1)
ppr=lag(ppr,1)
dn=lag(dn,1)

out.sbt=lm(resst0ta^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resst0ta^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resst0ta^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resst0ta^2 ~ dn + nnr + ppr)
summary(out.jsbt)

out.sbt=lm(resstgta^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resstgta^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resstgta^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resstgta^2 ~ dn + nnr + ppr)
summary(out.jsbt)

out.sbt=lm(resstata^2 ~ dn)
summary(out.sbt)
out.nsbt=lm(resstata^2 ~ nnr)
summary(out.nsbt)
out.psbt=lm(resstata^2 ~ ppr)
summary(out.psbt)
out.jsbt=lm(resstata^2 ~ dn + nnr + ppr)
summary(out.jsbt)

## ----fig.align='center'-------------------------------------------------------
ICall<-cbind(infocriteria(fit0n),infocriteria(fitgn),infocriteria(fitan),
infocriteria(fit0ta),infocriteria(fitgta),infocriteria(fitata))
colnames(ICall)<-c("GARCH N","GJR N","APARCH N","GARCH TA","GJR TA","APARCH TA")
ICall_df <- data.frame(ICall)
print(ICall_df)

## -----------------------------------------------------------------------------
p_EEM = log(window(EEMd, end="2023-03-31")) #log-prices
p_EEM = as.vector(p_EEM$EEM.Adjusted[,1])
T = length(p_EEM)
r_EEM=100*(p_EEM[2:T]-p_EEM[1:T-1]) #log percentage returns
y = r_EEM - mean(r_EEM)

# GARCH(1,1) Normal
for0n = ugarchroll(spec0n,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
for0nd = as.data.frame(for0n)
s0nf = for0nd$Sigma

# GARCH(1,1) Skewed-Student t
for0ta = ugarchroll(spec0ta,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
for0tad = as.data.frame(for0ta)
s0taf = for0tad$Sigma

# GJRGARCH(1,1) Normal
forgn = ugarchroll(specgn,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forgnd = as.data.frame(forgn)
sgnf = forgnd$Sigma

# GJRGARCH(1,1) Skewed-Student t
forgta = ugarchroll(specgta,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forgtad = as.data.frame(forgta)
sgtaf = forgtad$Sigma

# APARCH(1,1) Normal
foran = ugarchroll(specan,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
forand = as.data.frame(foran)
sanf = forand$Sigma

# APARCH(1,1) Skewed-Student t
forata = ugarchroll(specata,data=y,forecast.length=252,refit.every=5,window.size=1000, refit.window="moving")
foratad = as.data.frame(forata)
sataf = foratad$Sigma

## ----fig.align='center'-------------------------------------------------------
fall = cbind(s0nf,s0taf,sgnf,sgtaf,sanf,sataf)
par(mfrow=c(1,1))
matplot(fall, type = c("b"),pch=1,col = 1:6)
f0 = cbind(s0nf,s0taf)
fg = cbind(sgnf,sgtaf)
fa = cbind(sanf,sataf)
matplot(f0, type=('b'),pch=1,col=1:2)
matplot(fg, type=c('b'),pch=1,col=1:2)
matplot(fa, type=c('b'),pch=1,col=1:2)

## -----------------------------------------------------------------------------
rf = for0nd$Realized

lN = (s0nf^2 -rf^2)^2 # GARCH
lTA= (s0taf^2 -rf^2)^2 # GARCH skewed-t
lGN = (sgnf^2 -rf^2)^2 # GJRGARCH
lGTA = (sgtaf^2 -rf^2)^2 # GJRGARCH skewed-t
lAN = (sanf^2 -rf^2)^2 # GJRGARCH t
lATA = (sataf^2 -rf^2)^2 # APARCH

## -----------------------------------------------------------------------------
dNTA = lN-lTA 
dNGN = lN-lGN 
dNGTA = lN-lGTA 
dNAN = lN-lAN 
dNATA = lN-lATA 

dTAGN = lTA-lGN 
dTAGTA = lTA-lGTA 
dTAAN = lTA-lAN 
dTAATA = lTA-lATA 

dGNGTA = lGN-lGTA 
dGNAN = lGN-lAN 
dGNATA = lGN-lATA 

dGTAAN = lGTA-lAN 
dGTAATA = lGTA-lATA 

dANATA = lAN-lATA

## -----------------------------------------------------------------------------
m<-floor(0.75*((NROW(rf))^(1/3)))

x1<-as.vector(matrix(1,nrow=NROW(rf)))

VNTA = NeweyWest(lm(dNTA~x1-1),lag=m,prewhite=0)
VNGN = NeweyWest(lm(dNGN~x1-1),lag=m,prewhite=0)
VNGTA = NeweyWest(lm(dNGTA~x1-1),lag=m,prewhite=0)
VNAN = NeweyWest(lm(dNAN~x1-1),lag=m,prewhite=0)
VNATA = NeweyWest(lm(dNATA~x1-1),lag=m,prewhite=0)

VTAGN = NeweyWest(lm(dTAGN~x1-1),lag=m,prewhite=0)
VTAGTA = NeweyWest(lm(dTAGTA~x1-1),lag=m,prewhite=0)
VTAAN = NeweyWest(lm(dTAAN~x1-1),lag=m,prewhite=0)
VTAATA = NeweyWest(lm(dTAATA~x1-1),lag=m,prewhite=0)

VGNGTA = NeweyWest(lm(dGNGTA~x1-1),lag=m,prewhite=0)
VGNAN = NeweyWest(lm(dGNAN~x1-1),lag=m,prewhite=0)
VGNATA = NeweyWest(lm(dGNATA~x1-1),lag=m,prewhite=0)

VGTAAN = NeweyWest(lm(dGTAAN~x1-1),lag=m,prewhite=0)
VGTAATA = NeweyWest(lm(dGTAATA~x1-1),lag=m,prewhite=0)

VANATA = NeweyWest(lm(dANATA~x1-1),lag=m,prewhite=0)

DM<-matrix(0,nrow=6,ncol=6)
# Output -> robust variances of regressor

# loss = model in row minus model il column
colnames(DM)<-c("GARCHnorm","GARCH TA","GJRnorm","GJR ta","APARCHnorm","APARCH TA")
rownames(DM)<-c("GARCHnorm","GARCH TA","GJRnorm","GJR ta","APARCHnorm","APARCH TA")

DM[1, 2] <- mean(dNTA) / sqrt(VNTA)
DM[1, 3] <- mean(dNGN) / sqrt(VNGN)
DM[1, 4] <- mean(dNGTA) / sqrt(VNGTA)
DM[1, 5] <- mean(dNAN) / sqrt(VNAN)
DM[1, 6] <- mean(dNATA) / sqrt(VNATA)
DM[2, 3] <- mean(dTAGN) / sqrt(VTAGN)
DM[2, 4] <- mean(dTAGTA) / sqrt(VTAGTA)
DM[2, 5] <- mean(dTAAN) / sqrt(VTAAN)
DM[2, 6] <- mean(dTAATA) / sqrt(VTAATA)
DM[3, 4] <- mean(dGNGTA) / sqrt(VGNGTA)
DM[3, 5] <- mean(dGNAN) / sqrt(VGNAN)
DM[3, 6] <- mean(dGNATA) / sqrt(VGNATA)
DM[4, 5] <- mean(dGTAAN) / sqrt(VGTAAN)
DM[4, 6] <- mean(dGTAATA) / sqrt(VGTAATA)
DM[5, 6] <- mean(dANATA) / sqrt(VANATA)
DM

## -----------------------------------------------------------------------------
p_VT = log(window(VTd, end="2023-03-31")) #log-prices
p_VT = as.vector(p_VT$VT.Adjusted[,1])
T = length(p_VT)
r_VT=100*(p_VT[2:T]-p_VT[1:T-1]) #log percentage returns
yVT = r_VT - mean(r_VT)

specataVT <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)),
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                        distribution.model="sstd")
fitataVT <- ugarchfit(specataVT, r_VT)
resstataVT <- residuals(fitataVT, standardize=TRUE)
fitataVT@fit$robust.matcoef

p_EEM = log(window(EEMd, end="2023-03-31")) #log-prices
p_EEM = as.vector(p_EEM$EEM.Adjusted[,1])
T = length(p_EEM)
r_EEM=100*(p_EEM[2:T]-p_EEM[1:T-1]) #log percentage returns
yEEM = r_EEM - mean(r_EEM)

specataEEM <- ugarchspec(variance.model = list(model="apARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         distribution.model="sstd")
fitataEEM <- ugarchfit(specataEEM, r_EEM)
resstataEEM <- residuals(fitataEEM, standardize=TRUE)
fitataEEM@fit$robust.matcoef

## ----fig.align='center'-------------------------------------------------------
lista1=c()
for (i in 1:(min(length(yEEM),length(yVT))-62)){
  lista1[i]=cor(as.vector(residuals(fitataVT,standardize=T)[i:(i+61)]),
                as.vector(residuals(fitataEEM,standardize=T)[(i+1310):(i+1310+61)]))
}
plot(lista1)

lista2=c()
for (i in 1:(min(length(yEEM),length(yVT))-252)){
  lista2[i]=cor(as.vector(residuals(fitataVT,standardize=T)[i:(i+251)]),
                as.vector(residuals(fitataEEM,standardize=T)[(i+1310):(i+1310+251)]))
}
plot(lista2)

## ----fig.align='center'-------------------------------------------------------
par(mfrow = c(2,2))
for (i in seq(from = 1, to = 3653, by = 240)) {
  x <- as.vector(residuals(fitataVT, standardize = TRUE)[i:(i + 61)])
  y <- as.vector(residuals(fitataEEM, standardize = TRUE)[(i + 1310):(i + 1310 + 61)])
    start_date <- as.character(time(VTd)[i])
  end_date <- as.character(time(VTd)[i + 61])
    titolo <- paste("Correlation:", round(cor(x, y), 2), "\n", start_date, "to", end_date)
    plot(x, y, main = titolo, cex.main = 0.8)
}

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(1,1))
library(ggplot2)
residuals_VT <- residuals(fitataVT, standardize = TRUE)[1441:1502]
residuals_EEM <- residuals(fitataEEM, standardize = TRUE)[2751:2812]
data <- data.frame(residuals_VT, residuals_EEM)
ggplot(data) +
  geom_path(aes(x = seq_along(residuals_VT), y = residuals_VT), color = "blue") +
  geom_path(aes(x = seq_along(residuals_EEM), y = residuals_EEM), color = "red") +
  labs(y = "Residuals", title = "Correlation: 0.62 Period: 2014-03-18 --> 2014-06-13") +
  ylim(-3, 3) +
  theme_classic() +
  theme(axis.title.x = element_blank())

## ----fig.align='center'-------------------------------------------------------
par(mfrow=c(1,1))
library(ggplot2)
residuals_VT <- residuals(fitataVT, standardize = TRUE)[961:1023]
residuals_EEM <- residuals(fitataEEM, standardize = TRUE)[2271:2333]
data <- data.frame(residuals_VT, residuals_EEM)
ggplot(data) +
  geom_path(aes(x = seq_along(residuals_VT), y = residuals_VT), color = "blue") +
  geom_path(aes(x = seq_along(residuals_EEM), y = residuals_EEM), color = "red") +
  labs(y = "Residuals", title = "Correlation: 0.95 Period: 2011-05-05 to 2011-08-02") +
  ylim(-3, 3) +
  theme_classic() +
  theme(axis.title.x = element_blank())

