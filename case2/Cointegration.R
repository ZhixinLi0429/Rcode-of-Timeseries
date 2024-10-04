####################################################################
############          R code for Cointegration          ############
####################################################################

####################################################################
###########             Set-up for use                   ###########
####################################################################

rm(list = ls())   # Clear workspace

####################################################################
## Set input file directory: change to your own directory
####################################################################

## For Windows
setwd("C:\\users\\...") # replace the dots by your own path, use getwd() to see current path

## For Mac
## setwd("/users/...") # replace the dots by your own path, use getwd() to see current path

####################################################################
## Set output file directory: change to your own directory
## The output tables are created in this directory
####################################################################

## For Windows
output="C:\\users\\ymeersch\\Desktop"

## For Mac
## output="/users/ymeersch/Desktop"


## Install required packages (code is only executed if the packages are not yet installed)
if(!require(tseries)){install.packages("tseries")}
if(!require(urca)){install.packages("urca")}
if(!require(knitr)){install.packages("knitr")}
if(!require(tools)){install.packages("tools")}
if(!require(stargazer)){install.packages("stargazer")}
if(!require(xtable)){install.packages("xtable")}
if(!require(rmarkdown)){install.packages("rmarkdown")}
if(!require(dynlm)){install.packages("dynlm")}
if(!require(lmtest)){install.packages("lmtest")}

## Load required packages
library(urca)
library(stargazer)
library(xtable)
library(knitr)
library(tools)
library(rmarkdown)
library(dynlm)
library(lmtest)
source("TimeSeriesFunctions.R")

####################################################################
###########     Load data into the workspace             ###########
####################################################################
# Read from file
data         = read.table("C_Y_Z_database.csv", header = TRUE,sep = ";", skip = 1)
lnC          = ts(log(data[,"Consumption.expenditures.per.capita"]),start = c(1929, 1), frequency = 1)
lnY          = ts(log(data[,"Real.Disposable.Income.per.Capita"]),start = c(1929, 1), frequency = 1)
lnZ          = ts(log(data[,"Average.Wealth.per.capita"]),start = c(1929, 1), frequency = 1)
dlnC         = diff(lnC)
dlnY         = diff(lnY)
dlnZ         = diff(lnZ)

## Plot data in levels
pdf("lnC_lnY_lnZ_plot.pdf", onefile = T, paper = 'A4r',width = 0,height = 0)
par(mar = c(2,2,0,2)) 
plot(lnC, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA, las = 1, bty = "l", ylim = c(6,13))
lines(lnY, col = "green3", lwd = 2)
lines(lnZ, col = "blue3", lwd = 2, lty = 1)
abline(h = seq(6, 13, .5), col = "darkgrey", lty = 2)
legend("bottomright", legend = c("ln(Consumption)","ln(Income)","ln(Wealth)"), lty = c(1,1,1), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = FALSE,cex = 1.5,bty="n")
dev.off()


####################################################################
##################       Static Model             ##################
####################################################################
staticModel  = lm(lnC ~ lnY + lnZ)
stargazer(staticModel, type = "text")

staticResid  = ts(staticModel$residuals, start = 1929, frequency = 1)
staticFitted = ts(staticModel$fitted.values, start = 1929, frequency = 1)

pdf("static_fitted_resid.pdf", onefile = T, paper = 'A4r',width = 0,height = 0)
par(mar=c(5,2,2,3))
plot(lnC, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA,ylim = c(min(lnC,staticResid,staticFitted),max(lnC,staticResid,staticFitted)),
     las = 1, bty = "u")
lines(staticFitted, col = "green3", lwd = 2)
par(new = T)
plot(staticResid, type = "l", axes = F, xlab = NA, ylab = NA, col = "blue3", lwd = 2, lty = 3)
axis(side = 4,las = 1, at = seq(round(min(staticResid),1), round(max(staticResid),1), by = 0.1))
abline(h = 0, col = "darkgrey", lwd = 2)
par(new = T, xpd = NA, mar = c(par()$mar + c(-5,+3,0,0)))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab=NA,ylab=NA)
legend("bottom", legend = c("Actual","Fitted","Residuals (right axis)"), lty = c(1,1,3), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = TRUE,cex=0.9,bty="n")
dev.off()

## CRDW Test
DW_static = dwtest(staticModel)

## ADF test on the residuals from the static regression
selectARIMA(staticResid, 10, 0, 0)
ARmodel  = arima(staticResid, order = c(2,0,0), include.mean = FALSE)
LjungBox(ARmodel$residuals,20,2)

ADFresid = ur.df(staticResid, lags = 1, type = "none")
summary(ADFresid)

####################################################################
##################          ADL in levels         ##################
####################################################################
tmp = selectADL(lnC,lnY,lnZ,3,3,3)

ADLlevels = dynlm(lnC ~ L(lnC,1) + L(lnC,2) + lnY + L(lnY,1) + L(lnY,2) + lnZ + L(lnZ,1) + L(lnZ,2))
summary(ADLlevels)
plot(ADLlevels$residuals)
LjungBox(ADLlevels$residuals,20,2)

######################################################################################
#####  Error Correction Model with lagged residual from the static regression   ######
######################################################################################
#Resid_1    = lag(staticResid,-1)
#tmp = selectECM(lnC,lnY,lnZ,Resid_1,3,3,3)

ECMmodel  = dynlm(dlnC ~ L(dlnC,1) + dlnY + L(dlnY,1) + dlnZ + L(dlnZ,1) + L(staticResid,1) )
summary(ECMmodel)

ECMresid  = ECMmodel$residuals
ECMfitted = ECMmodel$fitted.values

#rmarkdown::render("ECM1.Rmd", output_file = "ECM1.tex")

pdf("ECM_resid.pdf", onefile = T, paper = 'A4r',width = 0,height = 0)
par(mar=c(5,2.5,2,3))
plot(dlnC, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA,ylim = c(min(dlnC,ECMresid,ECMfitted),max(dlnC,ECMresid,ECMfitted)),
     las = 1, bty = "u")
lines(ECMfitted, col = "green3", lwd = 2)
par(new = T)
plot(ECMresid, type = "l", axes = F, xlab = NA, ylab = NA, col = "blue3", lwd = 2, lty = 3)
axis(side = 4,las = 1, at = seq(round(min(ECMresid),2), round(max(ECMresid),2), by = 0.05))
abline(h = 0, col = "darkgrey", lwd = 2)
par(new = T, xpd = NA, mar = c(par()$mar + c(-5,+3,0,0)))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab=NA,ylab=NA)
legend("bottom", legend = c("Actual","Fitted","Residuals (right axis)"), lty = c(1,1,3), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = TRUE,cex=0.9,bty="n")
dev.off()

################################################################
#####  Error Correction Model with cointegrating vector   ######
################################################################
## Initial values for the parameter estimation. Those outside the error correction terms are initialized using the values from the first ECM model.
## Those within the error correction term use the coefficients from the static model as initialization.






alphainit  = coef(ECMmodel)[[7]]
theta1init = coef(ECMmodel)[[2]]
delta0init = coef(ECMmodel)[[3]]
delta1init = coef(ECMmodel)[[4]]
kappa0init = coef(ECMmodel)[[5]]
kappa1init = coef(ECMmodel)[[6]]

beta0init  = coef(staticModel)[[1]]
beta1init  = coef(staticModel)[[2]]
beta2init  = coef(staticModel)[[3]]

lnC_1    = lag(lnC,-1)
lnY_1    = lag(lnY,-1)
lnZ_1    = lag(lnZ,-1)
dlnC_1   = lag(dlnC,-1)
dlnY_1   = lag(dlnY,-1)
dlnZ_1   = lag(dlnZ,-1)
ECMdata  = data.frame(ts.union(dlnC, lnY, lnZ, lnC_1, lnY_1, lnZ_1, dlnC_1, dlnY, dlnY_1, dlnZ, dlnZ_1) )
formula  = dlnC ~ theta1*dlnC_1 + delta0*dlnY + delta1*dlnY_1 + kappa0*dlnZ + kappa1*dlnZ_1 + alpha*( lnC_1 - beta0 - beta1*lnY_1 - beta2*lnZ_1 )
ECMcoint = nls(formula, data = ECMdata, start = list(theta1 = theta1init, 
                                                     delta0 = delta0init, delta1 = delta1init, kappa0 = kappa0init, kappa1 = kappa1init, 
                                                     alpha = alphainit, beta0 = beta0init, beta1 = beta1init, beta2 = beta2init) )
summary(ECMcoint)
