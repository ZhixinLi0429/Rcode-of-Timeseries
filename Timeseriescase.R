setwd("C:\\Users\\PaulRomer\\Desktop\\Rworkfile\\Timeseries") 



if(!require(urca)){install.packages("urca")}
if(!require(tools)){install.packages("tools")}
if(!require(stargazer)){install.packages("stargazer")}
if(!require(dynlm)){install.packages("forecast")}
if(!require(Hmisc)){install.packages("Hmisc")}
if(!require(dynlm)){install.packages("dynlm")}
if(!require(lmtest)){install.packages("lmtest")}
if(!require(vars)){install.packages("vars")}
if(!require(VAR.etp)){install.packages("VAR.etp")}
if(!require(tseries)){install.packages("tseries")}
if(!require(knitr)){install.packages("knitr")}
if(!require(xtable)){install.packages("xtable")}
if(!require(rmarkdown)){install.packages("rmarkdown")}
if(!require(mise)){install.packages("mise")}


######################################################################################################


library(urca)            
library(stargazer)
library(tools)
library(Hmisc)
library(dynlm)
library(lmtest)  
library(vars)
library(VAR.etp)
library(xtable)
library(knitr)
library(forecast)
library(rmarkdown)
source("TimeSeriesFunctions.R")    


source("TimeSeriesFunctions.R") 

######################################################################################################
## Import dataset
######################################################################################################

data <- read.table("data.csv",header=TRUE,sep ="," ) 
data.ts =ts(data, start=c(1974, 1), end=c(2017, 12),frequency=12) 

Year       <- data$Year
Month      <- data$Month
OIL_P      <- data$OIL_P
OIL_PROD   <- data$OIL_PROD
OIL_STOCKS <- data$OIL_STOCKS
WORLD_IP   <- data$WORLD_IP
US_IP      <- data$US_IP
US_CPI     <- data$US_CPI
Surprises  <- data$Surprises
View(data)




US_IP =data.ts[,"US_IP"] 
US_IP_diff =diff(US_IP, differences=1)




########################
######Case 1############
########################


#1. plot US_IP in levels and US_IP in differences

pdf("US_IP.pdf", onefile=T, paper="A4r",width = 0,height = 0) 


par(mar=c(2,2,0,2)) 
plot(US_IP, ylab=NULL,xlab=NULL,lwd=2,col="darkblue",bty="l",col.axis="black")
abline(h=c(100,200,300,400), col="darkgrey", lty=2) 
dev.off() 


pdf("US_IP_diff.pdf", onefile=T, paper='A4r',width=0,height=0)
par(mar=c(2,2,0,2))
plot(US_IP_diff, ylab=NULL, xlab=NULL,lwd=2, col="darkblue",bty="l",col.axis="black")
abline(h=c(-2,-1,0,1,2,3),col="darkgrey",lty=2)
dev.off()

#From here onward, use the reduced period 1974:1-2015:12
Reduced_US_IP_diff <- ts(US_IP_diff,start=c(1974,1), 
                      end = c(2015,12), frequency = 12)






#2. Inspect the ACF and PACF for ∆ ln Yt

acf(Reduced_US_IP_diff)
pacf(Reduced_US_IP_diff)






#3. construct and estimate an ARMA model for ∆ ln Yt

AIC_Reduced_US_IP_diff <- aic_table(Reduced_US_IP_diff,6,6)
BIC_Reduced_US_IP_diff <- bic_table(Reduced_US_IP_diff,6,6)
print(AIC_Reduced_US_IP_diff)
print(BIC_Reduced_US_IP_diff)



ARMA5_5  <- arima(Reduced_US_IP_diff, order = c(5,0,5))
ARMA1_1  <- arima(Reduced_US_IP_diff, order = c(1,0,1))
AR3      <- arima(Reduced_US_IP_diff, order = c(3,0,0))

ARMA3_3  <- arima(Reduced_US_IP_diff, order = c(3,0,3))


LR_test_ARMA1_1_ARMA5_5 <- lrtest(ARMA1_1, ARMA5_5)
print(LR_test_ARMA1_1_ARMA5_5)

LR_test_ARMA3_3_ARMA5_5 <- lrtest(ARMA3_3, ARMA5_5)
print(LR_test_ARMA3_3_ARMA5_5)

LR_test_AR3_ARMA5_5 <- lrtest(AR3, ARMA5_5)
print(LR_test_AR3_ARMA5_5)







#4. Inspect the ACF and PACF for the residuals of your chosen ARMA model

residuals_ARMA5_5 <- residuals(ARMA5_5)
par(mfrow = c(2, 1)) 
acf(residuals_ARMA5_5, main = "ACF of ARMA5_5 Model Residuals")
pacf(residuals_ARMA5_5, main = "PACF of ARMA5_5 Model Residuals")

Qstatistic_arma5_5=LjungBox(ARMA5_5$residuals, 30, 10)
print(Qstatistic_arma5_5)






#5. Is your preferred ARMA model stable over time?


fit = Reduced_US_IP_diff - residuals_ARMA5_5

par(mar=c(5,5,2,5), mfrow = c(1,1))
plot(Reduced_US_IP_diff, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA,ylim = c(-6,3),
     las = 1, bty = "u")

lines(fit, col = "green3", lwd = 2)
abline(h = c(0), col = "darkgrey", lty = 2)
par(new = T)
plot(residuals_ARMA5_5, type = "l", axes = F, xlab = NA, 
     ylab = NA, col = "blue3", lwd = 2, lty = 3)
axis(side = 4)
par(new = T, xpd = NA, mar = c(par()$mar + c(-5,+3,0,0)))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab=NA,ylab=NA)
legend("bottom", legend = c("Actual","Fitted","Residuals (right axis)"), lty = c(1,1,3), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = TRUE,cex=0.9,bty="n")




se_residuals_ARMA5_5 <- residuals(ARMA5_5)-mean(residuals(ARMA5_5))
sigma_ARMA5_5 <- sqrt(ARMA5_5$sigma2)
cs_ARMA5_5 <- cumsum(se_residuals_ARMA5_5) / (sqrt(length(Reduced_US_IP_diff))*sigma_ARMA5_5)
cs_ARMA5_5.ts = ts(cs_ARMA5_5,start = c(1974,1), end = c(2015,12), frequency = 12)

par(mar=c(2,2,0,2))
plot(cs_ARMA5_5.ts, ylab = NULL, xlab = NULL,lwd = 3, col="darkblue", 
     bty = 'l', col.axis = "black", type = "l", 
     ylim = c(min(cs_ARMA5_5, -1.5), max(cs_ARMA5_5,1.5)))

abline(h = c(-1.22,1.22), col = "darkgrey", lty = 2,lwd = 2)
abline(h = c(-1.36,1.36), col = "grey", lty = 2,lwd = 2)
abline(h = c(-1.63,1.63), col = "green", lty = 2,lwd = 2)
abline(h = c(0), col = "black", lwd = 2)





#6. Use your preferred ARMA model to construct out-of-sample forecasts

predict_arma5_5_2016 <- predict(ARMA5_5, n.ahead = 12)
pred_2016 <- predict_arma5_5_2016$pred
plot(Reduced_US_IP_diff, type = "l", ylab = "US_IP_diff", 
     xlab = "Year", col = "blue")
lines(pred_2016, type="l", col="red")



actual_2016 <- ts(US_IP_diff,start=c(2016,1), end=c(2016,12),
                  frequency=12)
accuracy(pred_2016, actual_2016)




predict_arma5_5_2017 <- predict(ARMA5_5, n.ahead = 24)
pred_2017 <- predict_arma5_5_2017$pred
plot(Reduced_US_IP_diff, type = "l", ylab = "US_IP_diff", 
     xlab = "Year", col = "black")
lines(pred_2017, type="l", col="red")



actual_2017 <- ts(US_IP_diff,start=c(2017,1), end=c(2017,12),
                  frequency=12)
accuracy(pred_2017, actual_2017)




#7. Test whether ∆ lnYt is stationary

sample_size <- c(528)
n <-sample_size^(1/3)
print(n)


AR_AIC_US_IP_diff <- aic_table(US_IP_diff,10,0)
AR_BIC_US_IP_diff <- bic_table(US_IP_diff,10,0)
print(AR_AIC_US_IP_diff)
print(AR_BIC_US_IP_diff)

AR3 <- arima(US_IP_diff, order = c(3,0,0))
Qstatistic_ar3=LjungBox(AR3$residuals, 30, 3)
print(Qstatistic_ar3)

adf2_trend = ur.df (US_IP_diff, lags= 2, type= "trend")
summary(adf2_trend)



#8. Test whether lnYt is stationary

AR4_level      <- arima(US_IP, order = c(4,0,0))
Qstatistic_ar4_level =LjungBox(AR4_level$residuals, 30, 4)
print(Qstatistic_ar4_level )


level_adf3_trend = ur.df (US_IP, lags= 3, type= "trend")
summary(level_adf3_trend)

#US_IP_difflag1 <-lag(US_IP_diff, -1)
#US_IP_difflag2 <-lag(US_IP_diff, -2)
#US_IP_difflag3 <-lag(US_IP_diff, -3)
#trend <- seq(1, length(US_IP_diff))
#IfFtestsignificant <- lm(US_IP_diff~trend+US_IP_difflag1+US_IP_difflag2+US_IP_difflag3 )
#summary(IfFtestsignificant) it can also work if I need to remove the lnyt-1项 and do the regression.


level_adf3_drift = ur.df (US_IP, lags= 3, type= "drift")
summary(level_adf3_drift)

#US_IP_difflag1 <-lag(US_IP_diff, -1)
#US_IP_difflag2 <-lag(US_IP_diff, -2)
#US_IP_difflag3 <-lag(US_IP_diff, -3)
#IfFtestsignificant <- lm(US_IP_diff~US_IP_difflag1+US_IP_difflag2+US_IP_difflag3)
#summary(IfFtestsignificant) it can also work if I need to remove the lnyt-1项 and do the regression.


level_adf3_none = ur.df (US_IP, lags= 3, type= "none")
summary(level_adf3_none)
# here tau1 is 2.2615，critical value is -1.95，so do not reject the null hypothesis





####################################################
#################### PART 2 ########################
####################################################


#1. Construct and estimate an ADL model for the growth rate of U.S. industrial production ∆ ln Yt (use the full sample period)


Yt <- ts(US_IP, start=c(1974,1) ,frequency = 12)
diff_Yt <- diff(Yt, 1)

OIL_P <- ts(OIL_P, start=c(1974,1) ,frequency = 12)
diff_OIL_P <- diff(OIL_P, 1)

OIL_PROD <- ts(OIL_PROD, start=c(1974,1) ,frequency = 12)
diff_OIL_PROD <- diff(OIL_PROD, 1)

OIL_STOCKS <- ts(OIL_STOCKS, start=c(1974,1) ,frequency = 12)
diff_OIL_STOCKS <- diff(OIL_STOCKS, 1)

WORLD_IP <- ts(WORLD_IP, start=c(1974,1) ,frequency = 12)
diff_WORLD_IP <- diff(WORLD_IP, 1)

US_CPI <- ts(US_CPI, start=c(1974,1) ,frequency = 12)
diff_US_CPI <- diff(US_CPI, 1)


#construct and estimate the ADL model

aic_diffADL <- numeric()
for (i in 1:15) {
  adl_model <- dynlm(diff_Yt ~ lag(diff_Yt, i) 
                     + diff_OIL_P + lag(diff_OIL_P, i)
                     + diff_OIL_PROD + lag(diff_OIL_PROD, i)
                     + diff_OIL_STOCKS + lag(diff_OIL_STOCKS, i)
                     + diff_WORLD_IP + lag(diff_WORLD_IP, i)
                     + diff_US_CPI+ lag(diff_US_CPI, i)
                     , data = data)
  aic_diffADL[i] <- AIC(adl_model)
}
best_AICorder <- which.min(aic_diffADL)
print(best_AICorder)


Bic_diffADL <- numeric()
for (i in 1:15) {
  adl_model <- dynlm(diff_Yt ~ lag(diff_Yt, i) 
                     + diff_OIL_P + lag(diff_OIL_P, i)
                     + diff_OIL_PROD + lag(diff_OIL_PROD, i)
                     + diff_OIL_STOCKS + lag(diff_OIL_STOCKS, i)
                     + diff_WORLD_IP + lag(diff_WORLD_IP, i)
                     + diff_US_CPI+ lag(diff_US_CPI, i)
                     , data = data)
  Bic_diffADL[i] <- BIC(adl_model)
}
best_BICorder <- which.min(Bic_diffADL)
print(best_BICorder)



ADL_9 = dynlm(diff_Yt ~ L(diff_Yt,1:9) 
              + diff_OIL_P + L(diff_OIL_P,1:9)
              + diff_OIL_PROD + L(diff_OIL_PROD,1:9)
              + diff_OIL_STOCKS +L(diff_OIL_STOCKS, 1:9)
              + diff_WORLD_IP+ L(diff_WORLD_IP,1:9)
              + diff_US_CPI + L(diff_US_CPI,1:9)
)
stargazer(ADL_9,type = "text")


# diagnostic checking of the residuals
acf(ADL_9$residuals)
pacf(ADL_9$residuals)

LjungBox(ADL_9$residuals,80,59)



# Diagnostic checking of the parameter stability

fit_v = diff_Yt - ADL_9$residuals
par(mar=c(5,5,2,5), mfrow = c(1,1))
plot(diff_Yt, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA,ylim = c(-6,3),
     las = 1, bty = "u")
lines(fit_v, col = "green3", lwd = 2)
abline(h = c(0), col = "darkgrey", lty = 2)
par(new = T)
plot(ADL_9$residuals, type = "l", axes = F, xlab = NA, ylab = NA, col = "blue3", lwd = 2, lty = 3)
axis(side = 4)
par(new = T, xpd = NA, mar = c(par()$mar + c(-5,+3,0,0)))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab=NA,ylab=NA)
legend("bottom", legend = c("Actual","Fitted","Residuals (right axis)"), lty = c(1,1,3), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = TRUE,cex=0.9,bty="n")


se_residuals_ADL9<- residuals(ADL_9)-mean(residuals(ADL_9))
sigma_ADL_9 <- sigma(ADL_9)
cs_ADL_9 <- cumsum(se_residuals_ADL9) / (sqrt(length(diff_Yt))*sigma_ADL_9)
cs_ADL_9.ts = ts(cs_ADL_9,start = c(1974,1), end = c(2017,12), frequency = 12)

par(mar=c(2,2,0,2))
plot(cs_ADL_9.ts, ylab = NULL, xlab = NULL,lwd = 3, col="darkblue", bty = 'l', col.axis = "black", type = "l", ylim = c(min(cs_ADL_9, -1.5), max(cs_ADL_9,1.5)))
abline(h = c(-1.22,1.22), col = "darkgrey", lty = 2,lwd = 2)
abline(h = c(-1.36,1.36), col = "grey", lty = 2,lwd = 2)
abline(h = c(-1.63,1.63), col = "green", lty = 2,lwd = 2)
abline(h = c(0), col = "black", lwd = 2)


#2. Estimate a static model relating U.S. industrial production ln Yt to the levels of OIL P, OIL PROD, OIL STPOCKS,WORLD IP and US CPI

static_model_inlevel <- lm(Yt ~ OIL_P + OIL_PROD + OIL_STOCKS
                         +  WORLD_IP + US_CPI)
stargazer(static_model_inlevel, type = "text")


static_model_inlevel_residuals <-static_model_inlevel$residuals

pacf(static_model_inlevel_residuals)


AR2_residuals_static = arima(static_model_inlevel_residuals,
                         order=c(2,0,0), include.mean = TRUE)
AR4_residuals_static = arima(static_model_inlevel_residuals,
                         order=c(4,0,0), include.mean = TRUE)
LR_test_residuals_AR2_AR4 <- lrtest(AR2_residuals_static, AR4_residuals_static)
print(LR_test_residuals_AR2_AR4)



LjungBox(AR2_residuals_static$residuals, 25, 2)



ADF1_cointegration = ur.df(static_model_inlevel_residuals, lags = 1, type = "none")
summary(ADF1_cointegration)

MacKinnon_CriticalValue = -3.7429 - 8.352/528 - 13.41/528
print(MacKinnon_CriticalValue)



#3. Construct and estimate an ECM using the estimated residuals from your static regression

staticResiduals  = ts(static_model_inlevel$residuals, start = 1974, frequency = 12)

ECMmodel  = dynlm(diff_Yt ~ L(diff_Yt,1:9) 
                  + diff_OIL_P  + L(diff_OIL_P,1:9) 
                  + diff_OIL_PROD   + L(diff_OIL_PROD,1:9) 
                  + diff_OIL_STOCKS + L(diff_OIL_STOCKS,1:9)
                  + diff_WORLD_IP   + L(diff_WORLD_IP,1:9)
                  + diff_US_CPI     + L(diff_US_CPI,1:9)
                  + L(staticResiduals,1) )

### general diagonostic test for ECM

ECMresid  = ts(ECMmodel$residuals, start = 1974, frequency = 12)
ECMfitted = ts(ECMmodel$fitted.values, start = 1974, frequency = 12)
pdf("ECM_fitted_resid.pdf", onefile = T, paper = 'A4r',width = 0,height = 0)
par(mar=c(5,2,2,3))
plot(Yt, type = "l",col="red3", lwd = 2, xlab = NA, ylab = NA,ylim = c(min(Yt,ECMresid,ECMfitted),max(Yt,ECMresid,ECMfitted)),
     las = 1, bty = "u")
lines(ECMfitted, col = "green3", lwd = 2)
par(new = T)
plot(ECMresid, type = "l", axes = F, xlab = NA, ylab = NA, col = "blue3", lwd = 2, lty = 3)
axis(side = 4,las = 1, at = seq(round(min(ECMresid),1), round(max(ECMresid),1), by = 0.1))
abline(h = 0, col = "darkgrey", lwd = 2)
par(new = T, xpd = NA, mar = c(par()$mar + c(-5,+3,0,0)))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab=NA,ylab=NA)
legend("bottom", legend = c("Actual","Fitted","Residuals (right axis)"), lty = c(1,1,3), lwd = c(2,2,2)
       ,col = c("red3","green3","blue3"), horiz = TRUE,cex=0.9,bty="n")
dev.off()


se_ECMresid <- residuals(ECMmodel)-mean(residuals(ECMmodel))
sigma_ECMmodel <- sigma(ECMmodel)
cs_ECMmodel <- cumsum(se_ECMresid) / (sqrt(length(diff_Yt))*sigma_ECMmodel)
cs_ECMmodel.ts = ts(cs_ECMmodel,start = c(1974,1), end = c(2017,12), frequency = 12)

par(mar=c(2,2,0,2))
plot(cs_ECMmodel.ts, ylab = NULL, xlab = NULL,lwd = 3, col="darkblue", bty = 'l', col.axis = "black", type = "l", ylim = c(min(cs_ADL_9, -1.5), max(cs_ADL_9,1.5)))
abline(h = c(-1.22,1.22), col = "darkgrey", lty = 2,lwd = 2)
abline(h = c(-1.36,1.36), col = "grey", lty = 2,lwd = 2)
abline(h = c(-1.63,1.63), col = "green", lty = 2,lwd = 2)
abline(h = c(0), col = "black", lwd = 2)


LjungBox(ECMresid, 40,12)




stargazer(ECMmodel, type = "text")




#4. Construct an ADL model for ln Yt but now with all explanatory variables in levels


ADL_10 = dynlm(Yt ~ L(Yt,1:10) 
               + OIL_P + L(OIL_P,1:10)
               + OIL_PROD + L(OIL_PROD,1:10)
               + OIL_STOCKS +L(OIL_STOCKS, 1:10)
               + WORLD_IP+ L(WORLD_IP,1:10)
               + US_CPI + L(US_CPI,1:10)
)
stargazer(ADL_10,type = "text")


LjungBox(ADL_10$residuals, 100,65)




#5. Construct and estimate an ECM in which you estimate the long-run relationship along with the short-run dynamics

staticResiduals  = ts(static_model_inlevel$residuals, start = 1974, frequency = 12)

ECMmodel2  = dynlm(diff_Yt ~ L(diff_Yt,1) 
                  + diff_OIL_P  + L(diff_OIL_P,1) 
                  + diff_OIL_PROD   + L(diff_OIL_PROD,1) 
                  + diff_OIL_STOCKS + L(diff_OIL_STOCKS,1)
                  + diff_WORLD_IP   + L(diff_WORLD_IP,1)
                  + diff_US_CPI     + L(diff_US_CPI,1)
                  + L(staticResiduals,1) )
summary(ECMmodel2)



alphainit  = coef(ECMmodel2)[[13]]

theta1init = coef(ECMmodel2)[[2]]

delta0init = coef(ECMmodel2)[[3]]
delta1init = coef(ECMmodel2)[[4]]


kappa0init = coef(ECMmodel2)[[5]]
kappa1init = coef(ECMmodel2)[[6]]

Oilstock0init = coef(ECMmodel2)[[7]]
Oilstock1init = coef(ECMmodel2)[[8]]

worldip0init = coef(ECMmodel2)[[9]]
worldip1init = coef(ECMmodel2)[[10]]


uscpi0init = coef(ECMmodel2)[[11]]
uscpi1init = coef(ECMmodel2)[[12]]

beta0init  = coef(static_model_inlevel)[[1]]
beta1init  = coef(static_model_inlevel)[[2]]
beta2init  = coef(static_model_inlevel)[[3]]
beta3init  = coef(static_model_inlevel)[[4]]
beta4init  = coef(static_model_inlevel)[[5]]
beta5init  = coef(static_model_inlevel)[[6]]


Yt_1        = lag(Yt,-1)
OIL_P_1    = lag(OIL_P,-1)
OIL_PROD_1    = lag(OIL_PROD,-1)
OIL_STOCKS_1        = lag(OIL_STOCKS,-1)
WORLD_IP_1    = lag(WORLD_IP,-1)
US_CPI_1    = lag(US_CPI,-1)

diff_Yt_1   = lag(diff_Yt,-1)
diff_OIL_P_1   = lag(diff_OIL_P,-1)
diff_OIL_PROD_1   = lag(diff_OIL_PROD,-1)
diff_OIL_STOCKS_1   = lag(diff_OIL_STOCKS,-1)
diff_WORLD_IP_1   = lag(diff_WORLD_IP,-1)
diff_US_CPI_1   = lag(diff_US_CPI,-1)



ECMdata  = data.frame(ts.union(Yt, OIL_P, OIL_PROD, OIL_STOCKS, WORLD_IP, US_CPI,
                               Yt_1, OIL_P_1, OIL_PROD_1, OIL_STOCKS_1, WORLD_IP_1, US_CPI_1,
                               diff_Yt, diff_OIL_P, diff_OIL_PROD, diff_OIL_STOCKS, diff_WORLD_IP, diff_US_CPI,
                               diff_Yt_1, diff_OIL_P_1, diff_OIL_PROD_1, diff_OIL_STOCKS_1, diff_WORLD_IP_1, diff_US_CPI_1) )


formula  = diff_Yt ~ theta1*diff_Yt_1 + delta0*diff_OIL_P + delta1*diff_OIL_P_1 +
  kappa0*diff_OIL_PROD + kappa1*diff_OIL_PROD_1 + 
  alpha*( Yt_1 - beta0 - beta1*OIL_P_1 - beta2*OIL_PROD_1 )

ECMcoint = nls(formula, data = ECMdata, start = list(theta1 = theta1init, 
                                                     delta0 = delta0init, delta1 = delta1init, kappa0 = kappa0init, kappa1 = kappa1init,
                                                     alpha = alphainit, 
                                                     beta0 = beta0init, beta1 = beta1init, beta2 = beta2init))
summary(ECMcoint)




# the approach for ECM(9), full variables,there is singular matrix, very likely because of the multicolinearity



# Yt_lags <- lapply(1:9, function(i) lag(Yt, -i))
# OIL_P_lags <- lapply(1:9, function(i) lag(OIL_P, -i))
# OIL_PROD_lags <- lapply(1:9, function(i) lag(OIL_PROD, -i))
# OIL_STOCKS_lags <- lapply(1:9, function(i) lag(OIL_STOCKS, -i))
# WORLD_IP_lags <- lapply(1:9, function(i) lag(WORLD_IP, -i))
# US_CPI_lags <- lapply(1:9, function(i) lag(US_CPI, -i))


# ECMdata <- data.frame(Yt = Yt, Y_lags, OIL_P_lags, OIL_PROD_lags, OIL_STOCKS_lags, WORLD_IP_lags, US_CPI_lags)

# names(ECMdata) <- c("Yt", paste0("Y_lags", 1:9), paste0("OIL_P_lags", 1:9), 
                    # paste0("OIL_PROD_lags", 1:9), paste0("OIL_STOCKS_lags", 1:9), 
                    # paste0("WORLD_IP_lags", 1:9), paste0("US_CPI_lags", 1:9))


# formula <- as.formula(paste("Yt ~",paste(c(paste0("Y_lags", 1:9),
                                           # paste0("OIL_P_lags", 1:9),
                                           # paste0("OIL_PROD_lags", 1:9),
                                           # paste0("OIL_STOCKS_lags", 1:9),
                                           # paste0("WORLD_IP_lags", 1:9),
                                           # paste0("US_CPI_lags", 1:9)), collapse = "+")))

# initial_model <- ECMmodel
# initial_values <- coef(initial_model)

# ECMcointegration <- nls(formula, data=ECMdata, start = as.list(initial_values))
# summary(ECMcointegration)




# the approach for ECM(1), full variables,there is still singular gradient, i.e. rank deficiency, very likely because of the multicolinearity

#alphainit  = coef(ECMmodel2)[[13]]
#theta1init = coef(ECMmodel2)[[2]]
#delta0init = coef(ECMmodel2)[[3]]
#delta1init = coef(ECMmodel2)[[4]]
#kappa0init = coef(ECMmodel2)[[5]]
#kappa1init = coef(ECMmodel2)[[6]]
#Oilstock0init = coef(ECMmodel2)[[7]]
#Oilstock1init = coef(ECMmodel2)[[8]]
#worldip0init = coef(ECMmodel2)[[9]]
#worldip1init = coef(ECMmodel2)[[10]]
#uscpi0init = coef(ECMmodel2)[[11]]
#uscpi1init = coef(ECMmodel2)[[12]]
#beta0init  = coef(static_model_inlevel)[[1]]
#beta1init  = coef(static_model_inlevel)[[2]]
#beta2init  = coef(static_model_inlevel)[[3]]
#beta3init  = coef(static_model_inlevel)[[4]]
#beta4init  = coef(static_model_inlevel)[[5]]
#beta5init  = coef(static_model_inlevel)[[6]]


#Yt_1        = lag(Yt,-1)
#OIL_P_1    = lag(OIL_P,-1)
#OIL_PROD_1    = lag(OIL_PROD,-1)
#OIL_STOCKS_1        = lag(OIL_STOCKS,-1)
#WORLD_IP_1    = lag(WORLD_IP,-1)
#US_CPI_1    = lag(US_CPI,-1)

#diff_Yt_1   = lag(diff_Yt,-1)
#diff_OIL_P_1   = lag(diff_OIL_P,-1)
#diff_OIL_PROD_1   = lag(diff_OIL_PROD,-1)
#diff_OIL_STOCKS_1   = lag(diff_OIL_STOCKS,-1)
#diff_WORLD_IP_1   = lag(diff_WORLD_IP,-1)
#diff_US_CPI_1   = lag(diff_US_CPI,-1)

#ECMdata  = data.frame(ts.union(Yt, OIL_P, OIL_PROD, OIL_STOCKS, WORLD_IP, US_CPI,
#Yt_1, OIL_P_1, OIL_PROD_1, OIL_STOCKS_1, WORLD_IP_1, US_CPI_1,
#diff_Yt, diff_OIL_P, diff_OIL_PROD, diff_OIL_STOCKS, diff_WORLD_IP, diff_US_CPI,
#diff_Yt_1, diff_OIL_P_1, diff_OIL_PROD_1, diff_OIL_STOCKS_1, diff_WORLD_IP_1, diff_US_CPI_1) )


#formula  = diff_Yt ~ theta1*diff_Yt_1 + delta0*diff_OIL_P + delta1*diff_OIL_P_1 +
#kappa0*diff_OIL_PROD + kappa1*diff_OIL_PROD_1 + Oilstock0*diff_OIL_PROD + Oilstock1*diff_OIL_PROD_1+
#worldip0*diff_WORLD_IP + worldip1*diff_WORLD_IP_1 + uscpi0*diff_US_CPI + uscpi1*diff_US_CPI_1+ 
#alpha*( Yt_1 - beta0 - beta1*OIL_P_1 - beta2*OIL_PROD_1 - beta3*OIL_STOCKS_1 - beta4*WORLD_IP_1 - beta5*US_CPI_1 )

#ECMcoint = nls(formula, data = ECMdata, start = list(theta1 = theta1init, 
#delta0 = delta0init, delta1 = delta1init, kappa0 = kappa0init, kappa1 = kappa1init,
#Oilstock0 = Oilstock0init, Oilstock1 = Oilstock1init, 
#worldip0 = worldip0init, worldip1 = worldip1init, uscpi0 = uscpi0init, uscpi1 = uscpi1init,
#alpha = alphainit, 
#beta0 = beta0init, beta1 = beta1init, beta2 = beta2init, beta3 = beta3init,
#beta4 = beta4init, beta5 = beta5init) )
#summary(ECMcoint)










######################################################################################################
################################################ CASE 3 ##############################################
######################################################################################################


if(!require(vars)){install.packages("vars")}
library(vars)
if(!require(tidyr)){install.packages("tidyr")}
library(tidyr)
if(!require(dplyr)){install.packages("dplyr")}
library(dplyr)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
if(!require(devtools)){install.packages("devtools")}
library(devtools)
if(!require(tidyverse)){install.packages("tidyverse")}
library(tidyverse)
remotes::install_github("angusmoore/varexternalinstrument")
library(varexternalinstrument)


library(VAR.etp)
library(varexternalinstrument)
library(varexternal)
library(vars)




#1. Estimate a reduced-form VAR including OIL P, OIL PROD,OIL STOCKS, WORLD IP, US IP and US CPI



VARselect(data_var, lag.max=20, type = "const")



data_var=cbind(Yt,OIL_P,OIL_PROD,OIL_STOCKS,WORLD_IP,US_CPI)
Reduced_VAR12 = VAR(data_var, p = 12, type = "const")
summary(Reduced_VAR12)


plot(Reduced_VAR12)
serialcorrelation(VARmodel = Reduced_VAR12,nlag=25)


var12.stab=stability(Reduced_VAR12,type = "OLS-CUSUM")
plot(var12.stab, alpha=0.05)





#2. Perform a Granger Causality analysis


#Joint-test


causality(Reduced_VAR12, cause = c("OIL_P", "OIL_PROD","OIL_STOCKS",
                                   "WORLD_IP","US_CPI"))$Granger

causality(Reduced_VAR12, cause = c("Yt","OIL_PROD","OIL_STOCKS",
                                   "WORLD_IP","US_CPI"))$Granger

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_STOCKS",
                                   "WORLD_IP","US_CPI"))$Granger

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD",
                                   "WORLD_IP","US_CPI"))$Granger

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD",
                                   "OIL_STOCKS","US_CPI"))$Granger

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD",
                                   "OIL_STOCKS","WORLD_IP"))$Granger


#univariate-test

causality(Reduced_VAR12, cause = c("Yt"))$Granger

causality(Reduced_VAR12, cause = c("OIL_P"))$Granger

causality(Reduced_VAR12, cause = c("OIL_PROD"))$Granger

causality(Reduced_VAR12, cause = c("OIL_STOCKS"))$Granger

causality(Reduced_VAR12, cause = c("WORLD_IP"))$Granger

causality(Reduced_VAR12, cause = c("US_CPI"))$Granger






#  3. Use a Cholesky decomposition to identify the SVAR


##correlation matrix of the residuals


VARsum   = summary(Reduced_VAR12)
residCorr = VARsum$corres
stargazer(residCorr, type = "text")  

causality(Reduced_VAR12, cause = c("OIL_P", "OIL_PROD","OIL_STOCKS", "WORLD_IP","US_CPI"))$Instant

causality(Reduced_VAR12, cause = c("Yt","OIL_PROD","OIL_STOCKS", "WORLD_IP","US_CPI"))$Instant

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_STOCKS", "WORLD_IP","US_CPI"))$Instant

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD", "WORLD_IP","US_CPI"))$Instant

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD", "OIL_STOCKS","US_CPI"))$Instant

causality(Reduced_VAR12, cause = c("Yt","OIL_P","OIL_PROD", "OIL_STOCKS","WORLD_IP"))$Instant




VARord = cbind(OIL_PROD, OIL_P, OIL_STOCKS, US_CPI, Yt, WORLD_IP)

SVAR = VAR(VARord, p = 12, type = "const")

# 在第一列的OIL_PROD上引入10%的冲击

impact_index = 1  
impact_size = 0.10  

# 创建一个矩阵表示冲击

impact_matrix = matrix(0, ncol = length(VARord), nrow = length(VARord))
impact_matrix[impact_index, impact_index] = impact_size

# 应用冲击

SVAR$A = SVAR$A + impact_matrix

IRF_SVAR = irf(SVAR, n.ahead = 50, boot='TRUE',runs=100) 

plot(IRF_SVAR)



VD_SVAR<-fevd(SVAR, n.ahead = 50)

png("VD_SVAR.png", width = 800, height = 600)  # 根据需要设置宽度和高度

par(mar = c(5, 4, 4, 2) + 0.1)  # 适当调整边距参数

plot(VD_SVAR)

par(mar = c(5, 4, 4, 2))

dev.off()



#  4. Use the oil supply surprises variable as an external instrument for the real oil price to identify the structural VAR


#specify variables and instrument

new_data <- data.frame(Yt, OIL_P, OIL_PROD, OIL_STOCKS, WORLD_IP, US_CPI, Surprises)

y <- cbind(new_data$Yt,new_data$OIL_P,new_data$OIL_PROD, new_data$OIL_STOCKS,new_data$WORLD_IP,new_data$US_CPI)

colnames(y) <- c("Yt","OIL_P", "OIL_PROD", "OIL_STOCKS", "WORLD_IP", "US_CPI")

Z <- ts(Surprises)

names(Z) <- c("Surprises")


# estimate reduced var

GKvar <- VAR(y, p=12, type= "const")
shockcol <- externalinstrument(GKvar,Z, "OIL_P")
shockcol


# IRF

ma_reprensentation <- Phi(GKvar, 50)
irfs <- apply (ma_reprensentation, 3, function (x) x %*% shockcol)
irfs <- as.data.frame(t(irfs))
colnames(irfs) <- names(shockcol)
irfs <- mutate(irfs, horizon = 0:50)
irfs <- gather(irfs, key = variable, value = response, -horizon)
ggplot(irfs, aes(x=horizon, y=response, group= variable, color= variable))+geom_line()




irf_OIL_P <- ggplot((subset(irfs, variable %in% c("OIL_P"))), 
                    aes(x = horizon, y = response, group = variable, 
                        color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))
irf_OIL_PROD <- ggplot((subset(irfs, variable %in% c("OIL_PROD"))), 
                       aes(x = horizon, y = response, group = variable, 
                           color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))
irf_OIL_STOCKS <- ggplot((subset(irfs, variable %in% c("OIL_STOCKS"))), 
                         aes(x = horizon, y = response, group = variable, 
                             color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))
irf_WORLD_IP <- ggplot((subset(irfs, variable %in% c("WORLD_IP"))), 
                       aes(x = horizon, y = response, group = variable, 
                           color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))
irf_Yt <- ggplot((subset(irfs, variable %in% c("Yt"))), 
                 aes(x = horizon, y = response, group = variable, 
                     color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))
irf_US_CPI <- ggplot((subset(irfs, variable %in% c("US_CPI"))), 
                     aes(x = horizon, y = response, group = variable, 
                         color = variable)) + geom_line()+scale_x_continuous(breaks=seq(0,50, 5))



install.packages("cowplot")
library(cowplot)



irfs_plots <- plot_grid(irf_OIL_P, irf_OIL_PROD, irf_OIL_STOCKS, 
                        irf_WORLD_IP, irf_Yt, irf_US_CPI, 
                        ncol = 2, labels = "AUTO")
irfs_plots





#5. Use local projections as an alternative approach to obtain the response of each of the variables to an oil supply surprise



install.packages("lpirfs")
library(lpirfs)


##Local projection 


y <- cbind(new_data$Yt,new_data$OIL_P,new_data$OIL_PROD, new_data$OIL_STOCKS,new_data$WORLD_IP,new_data$US_CPI)

colnames(y) <- c("Yt","OIL_P", "OIL_PROD", "OIL_STOCKS", "WORLD_IP", "US_CPI")

y_lp <-as.data.frame(y)

LP = lp_lin(y_lp, lags_endog_lin = 12, shock_type = 1, trend = 0, confint = 1.96, hor = 50)
plot(LP) 


##lp with IV by 2SLS estimation

shock <- as.data.frame(new_data$OIL_P)

instrum <-as.data.frame(new_data$Surprises)

LP_iv = lp_lin_iv(y_lp, lags_endog_lin = 12, shock = shock, 
                  instrum = instrum, use_twosls = TRUE, trend = 0,
                  confint = 1.96, hor= 50)

plot(LP_iv)



