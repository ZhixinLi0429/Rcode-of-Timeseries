####################################################################
##################    Box-Pierce Q-statistic      ##################
####################################################################
LjungBox = function (object,nlag, fitdf) 
{
  if(!is.ts(object)) {
    rs     = object$model
    stdres = rs / sqrt(object$sigma2)
    #fitdf  = sum(object$arma[1:2])
    pval   = numeric(nlag - fitdf - 1)
    tstat  = numeric(nlag - fitdf - 1)
    for (i in ((fitdf + 1):nlag)) {
      pval[i-fitdf]  = round(Box.test(rs, i, type = "Ljung-Box",fitdf)$p.value, 4)
      if(pval[i-fitdf] < 0.0001) 
      {
        pval[i-fitdf] = "<0.0001"
      }
      tstat[i-fitdf] = round(Box.test(rs, i, type = "Ljung-Box",fitdf)$statistic, 4)
    }
  } else {
    pval   = numeric(nlag)
    tstat  = numeric(nlag)
    for (i in 1:nlag) {
      if(fitdf >= i) {
        pval[i] = NA
      } else {
        pval[i]  = round(Box.test(object, i, type = "Ljung-Box")$p.value, 4)
        if(pval[i] < 0.0001) {
          pval[i] = "<0.0001"
        }
      }
      tstat[i] = round(Box.test(object, i, type = "Ljung-Box")$statistic, 4) 
    }
  }
  pval           = data.frame(seq(1,nlag,1), tstat, pval)
  names(pval)    = c("lag","Ljung-Box","p-value")
  return(pval)
}

####################################################################
##################         IC tables              ##################
####################################################################
# Estimation method: CSS-ML
aic_table <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      tmp = arima(data,order=c(p,0,q), method = "CSS-ML", include.mean = TRUE, transform.pars = FALSE)
      table[p+1,q+1] <- AIC(tmp)
    }
  }
  dimnames(table) <- list(paste("AR",0:P, "", sep=""),paste("MA",0:Q,sep=""))
  table
}

bic_table <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      tmp = arima(data,order=c(p,0,q), method = "CSS-ML", include.mean = TRUE, transform.pars = FALSE)
      bic = AIC(tmp,k = log(length(data)))
      table[p+1,q+1] <- bic
    }
  }
  dimnames(table) <- list(paste("AR",0:P, "", sep=""),paste("MA",0:Q,sep=""))
  table
}

# Estimation: CSS 
aic_table_css <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      tmp = arima(data,order=c(p,0,q), method = "CSS", include.mean = TRUE, transform.pars = FALSE, n.cond = P)
      rss = sum(residuals(tmp)^2)
      aic = (tmp$nobs-P)*log(rss)+2*length(tmp$coef)
      table[p+1,q+1] <- aic
    }
  }
  dimnames(table) <- list(paste("AR",0:P, "", sep=""),paste("MA",0:Q,sep=""))
  table
}

bic_table_css <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      tmp = arima(data,order=c(p,0,q), method = "CSS", include.mean = TRUE, transform.pars = FALSE, n.cond = P)
      rss = sum(residuals(tmp)^2)
      bic = (tmp$nobs-P)*log(rss)+length(tmp$coef)*log(tmp$nobs-P)
      table[p+1,q+1] <- bic
    }
  }
  dimnames(table) <- list(paste("AR",0:P, "", sep=""),paste("MA",0:Q,sep=""))
  table
}

####################################################################
##################   AIC, BIC and RSS for NLS     ##################
####################################################################
####################################################################
IC = function (model) {
  adjT       = model$nobs
  k          = dim(model$var.coef)[1]
  RSS        = t(model$residuals)%*%model$residuals
  ll         = - (adjT/2)*( 1 + log(2*pi) + log(RSS / adjT) )
  BIC        = - 2*(ll/adjT) + (k/adjT)*log(adjT)
  AIC        = - 2*(ll/adjT) + (k/adjT)*2
  out        = list(round(AIC, 4), round(BIC, 4), round(ll, 4), round(RSS, 4) )
  names(out) = c("AIC", "BIC", "LogLik", "RSS")
  #stargazer(data.frame(out), type = "text", summary = FALSE)
  return(out)
}

## Automatically selects the best ADL model (in levels)

selectADL = function(var1, var2, var3, var4, var5, var6, k1, k2, k3, k4, k5, k6)
{
  # Table to capture AIC & BIC 
  AIC.table =  array(dim = c(k1+1, k2+1, k3+1, k4+1, k5+1, k6+1))
  BIC.table =  array(dim = c(k1+1, k2+1, k3+1, k4+1, k5+1, k6+1))
  if (max(k1, k2, k3, k4, k5, k6) > 10) {
    month_newvar = start(var1)[2] + max(k1, k2, k3, k4, k5, k6) - 12
 } else{
    month_newvar = start(var1)[2] + max(k1, k2, k3, k4, k5, k6)
  }
  
  if(max(k1, k2, k3, k4, k5, k6) < 11){
    year_newvar = start(var1)[1]
    
  }else if(max(k1, k2, k3, k4, k5, k6) < 23){
    year_newvar = start(var1)[1] + 1
    }else{
      year_newvar = start(var1)[1] + 2
      }
  
  newvar = window(var1, start = c(year_newvar,month_newvar), end = end(var1), frequency = frequency(var1))
  
  # Load AIC & BIC in table 
  for(i in 0:k1){
    for(j in 0:k2){
      for(k in 0:k3){
        for(l in 0:k4){
          for(m in 0:k5){
            for(n in 0:k6){
              if(i == 0){
                AIC.table[i+1, j+1, k+1, l+1, m+1, n+1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var2, 0:j) + L(var3, 0:k) + L(var4, 0:l) + L(var5, 0:m) + L(var6, 0:n), start = c(year_newvar,month_newvar)))^2)) + 2*(6 + i + j + k + l + m + n), digits = 4)
                BIC.table[i+1, j+1, k+1, l+1, m+1, n+1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var2, 0:j) + L(var3, 0:k) + L(var4, 0:l) + L(var5, 0:m) + L(var6, 0:n), start = c(year_newvar,month_newvar)))^2)) + (6 + i + j + k + l + m + n)*log(length(newvar)), digits = 4)
              }
              else{
              AIC.table[i+1, j+1, k+1,l+1, m+1, n+1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var1, 1:i) + L(var2, 0:j) + L(var3, 0:k) + L(var4, 0:l) + L(var5, 0:m) + L(var6, 0:n), start = c(year_newvar,month_newvar)))^2)) + 2*(6 + i + j + k + l + m + n), digits = 4)
              BIC.table[i+1, j+1, k+1,l+1, m+1, n+1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var1, 1:i) + L(var2, 0:j) + L(var3, 0:k) + L(var4, 0:l) + L(var5, 0:m) + L(var6, 0:n), start = c(year_newvar,month_newvar)))^2)) + (6 + i + j + k + l + m + n) * log(length(newvar)), digits = 4)
              }
            }
          } 
        }
      }
    }
  } 
  
  # 10% quantile for the AIC table 
  AICQ10 = quantile(AIC.table, 0.1, na.rm = TRUE)
  BICQ10 = quantile(BIC.table, 0.1, na.rm = TRUE)
  
  # Show the results for the 10% best 
  for (i in 0:k1){
    for(j in 0:k2){
      for(k in 0:k3){
        for(l in 0:k4){
          for(m in 0:k5){
            for(n in 0:k6){
              if(AIC.table[i+1, j+1, k+1,l+1, m+1, n+1] < AICQ10 || BIC.table[i+1, j+1, k+1,l+1, m+1, n+1] < BICQ10) {
                cat("Model (k1,k2,k3,k4,k5,k6) : ",i ,j,k,l,m,n, ":  nb param: ", 6+i+j+k+l+m+n, "    AIC:", AIC.table[i+1, j+1, k+1,l+1, m+1, n+1], "    BIC:", BIC.table[i+1, j+1, k+1,l+1, m+1, n+1],  "\n")
              }
            }
          }
        }
      }
    }
  }
}


## Automatically selects the best ECM model 

selectECM = function(var1, var2, var3, staticResid, k1, k2, k3)
{
  # Table to capture AIC & BIC 
  AIC.table =  array(dim = c(k1+1, k2+1, k3+1))
  BIC.table =  array(dim = c(k1+1, k2+1, k3+1))
  newstart = start(var1) + max(k1, k2, k3)
  newvar = window(var1, start = newstart[1], end = end(var1), frequency = frequency(var1))
  
  # Load AIC & BIC in table 
  for(i in 0:k1){
    for(j in 0:k2){
      for(k in 0:k3){
        if(i == 0){
          AIC.table[i + 1, j + 1, k + 1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var2, 0:j) + L(var3, 0:k) + L(staticResid, 1), start = newstart[1]))^2)) + 2*(4 + i + j + k), digits = 4)
          BIC.table[i + 1, j + 1, k + 1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var2, 0:j) + L(var3, 0:k) + L(staticResid, 1), start = newstart[1]))^2)) + (4 + i + j + k)*log(length(newvar)), digits = 4)
        }
        else{
          AIC.table[i + 1, j + 1, k + 1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var1, 1:i) + L(var2, 0:j) + L(var3, 0:k) + L(staticResid, 1), start = newstart[1]))^2)) + 2*(4 + i + j + k), digits = 4)
          BIC.table[i + 1, j + 1, k + 1] = round(length(newvar)*log(sum(residuals(dynlm(var1 ~ L(var1, 1:i) + L(var2, 0:j) + L(var3, 0:k) + L(staticResid, 1), start = newstart[1]))^2)) + (4 + i + j + k) * log(length(newvar)), digits = 4)
        }
      }
    }
  } 
  # 10% quantile for the AIC table 
  AICQ10 = quantile(AIC.table, 0.1, na.rm = TRUE)
  BICQ10 = quantile(BIC.table, 0.1, na.rm = TRUE)
  
  # Show the results for the 10% best 
  for (i in 0:k1){
    for(j in 0:k2){
      for(k in 0:k3){
        if(AIC.table[i + 1, j + 1, k + 1] < AICQ10 || BIC.table[i + 1, j + 1, k + 1 ] < BICQ10) {
          cat("Model (k1,k2,k3) : ",i ,j,k, ":  nb param: ", 4+i+j+k, "    AIC:", AIC.table[i + 1,j + 1, k + 1], "    BIC:", BIC.table[i + 1,j + 1, k + 1],  "\n")
        }
      }
    }
  }
}


## Automatically selects the 10% best models in terms of SBIC and AIC
selectARIMA = function(data, p.max, d, q.max){
  # Table to capture AIC & BIC
  AIC.table =  array(dim = c(p.max+1,q.max+1))
  BIC.table =  array(dim = c(p.max+1,q.max+1))
  
  # Load AIC & BIC in table
  for (i in 0:p.max) {
    for (j in 0:q.max) {
      AIC.table[i + 1,j + 1] = IC(arima(data, order = c(i,d, j)))$AIC
      BIC.table[i + 1,j + 1] = IC(arima(data, order = c(i,d, j)))$BIC
    }
  }
  
  # 10%-quantile for the AIC table
  AICQ10 = quantile(AIC.table, 0.1, na.rm= TRUE)
  BICQ10 = quantile(BIC.table, 0.1, na.rm= TRUE)       
  
  # Show results for the 10% best
  for (i in 0:p.max) {
    for (j in 0:q.max) {         
      if (AIC.table[i + 1,j + 1] < AICQ10 || BIC.table[i + 1,j + 1] < BICQ10) { 
        cat("Model (p,d,q) : ",i ,d, j, ":  nb param: ", i+j, "    AIC:", AIC.table[i + 1,j + 1], "    BIC:", BIC.table[i + 1,j + 1],  "\n") } }
    
  }
}



####################################################################
##################  Check whether object is empty  #################
####################################################################
isEmpty <- function(x) {
  return(ifelse(length(x)==0 || all(x==0), TRUE, FALSE))
}

####################################################################
######     Impulse response functions for ARIMA-models     #########
####################################################################
IRFarma = function(object, length) {
  IRF    = rep(as.numeric(object$coef["intercept"]), length)
  IRF[5] = IRF[5] + 1 # shock
  if(!isEmpty(object$model$phi) && isEmpty(object$model$theta)) {
    AR   = object$model$phi
    P    = length(AR)
    
    for (i in 1:(length - 5)) {
      for (p in 1:P) {
        if(i >= p) {
          IRF[5 + i] = IRF[5 + i] + t(AR[p])*(IRF[5 + i - p] - IRF[4])
          #IRF[5 + i] = IRF[5 + i] + t(AR[p]^(i - p + 1))
          #IRF[5 + i] = IRF[5 + i] + t(AR[p])*(IRF[4 + i] - IRF[4])
        }
      }
    }
    return(IRF)
  }
  if(isEmpty(object$model$phi) && !isEmpty(object$model$theta)) {
    MA   = object$model$theta
    Q    = length(MA)
    
    IRF[6:(6 + Q - 1)] = IRF[6:(6 + Q - 1)] + MA[1:Q]
    return(IRF)
  }
  if(!isEmpty(object$model$phi) && !isEmpty(object$model$theta)) {
    AR   = object$model$phi
    MA   = object$model$theta
    P    = length(AR)
    Q    = length(MA)
    
    IRF[6:(6 + Q - 1)] = IRF[6:(6 + Q - 1)] + MA[1:Q]
    
    for (i in 1:(length - 5)) {
      for (p in 1:P) {
        if(i >= p) {
          IRF[5 + i] = IRF[5 + i] + t(AR[p])*(IRF[5 + i - p] - IRF[4])
          #IRF[5 + i] = IRF[5 + i] + t(AR[p])*(IRF[5 + i - p] - IRF[5 + i])
          #IRF[5 + i] = IRF[5 + i] + t(AR[p]^(i - p + 1))
        }
      }
    }
    return(IRF)
  }
}

####################################################################
##################       Static forecast       #####################
####################################################################
staticforecast = function(data.ts, order, forecastLength){
  n         = length(data.ts)
  Pred      = rep(NA,n)
  startTime = time(data.ts)[1]
  freqTime  = frequency(time(data.ts))
  
  for(i in (n - forecastLength - 1):n){
    ts.part    = data.ts[1:(i-1)]
    tmp.model  = arima(ts.part, order = order)
    Pred[i]    = predict(tmp.model, n.ahead = 1)$pred[1]
  }
  tspred = ts(Pred, start = startTime, frequency = freqTime)
  
  return(tspred)
}


####################################################################
#######    Portmanteau, adjusted Portmanteau and BG-lm    ##########
####################################################################
serialcorrelation = function (VARmodel, nlag) 
{
  PTasympStat = numeric(nlag)
  PTadjTest   = numeric(nlag)
  LMstat      = numeric(nlag)
  
  PTasympPval = numeric(nlag)
  PTadjPval   = numeric(nlag)
  LMpval      = numeric(nlag)
  for (i in 1:nlag) {
    tempPTasymp    = serial.test(VARmodel, lags.pt = i, type = "PT.asymptotic")$serial
    tempPTadj      = serial.test(VARmodel, lags.pt = i, type = "PT.adjusted")$serial
    tempLM         = serial.test(VARmodel, lags.bg = i, type = "BG")$serial
    
    PTasympStat[i] = round(tempPTasymp$statistic, 4)
    PTadjTest[i]   = round(tempPTadj$statistic, 4)
    LMstat[i]      = round(tempLM$statistic, 4)
    
    PTasympPval[i] = round(tempPTasymp$p.value, 4)
    PTadjPval[i]   = round(tempPTadj$p.value, 4)
    LMpval[i]      = round(tempLM$p.value, 4)
    
    if(i<1+VARmodel$p) {PTasympPval[i] = "NaN"} #else {{if(PTasympPval[i] < 0.0001) PTasympPval[i] = "<0.0001"}}
    if(i<1+VARmodel$p) PTadjPval[i] = "NaN"
  }
  output           = data.frame(seq(1,nlag,1), PTasympStat, PTasympPval, PTadjTest, PTadjPval, LMstat, LMpval)
  names(output)    = c("lag","Portm. stat","Portm. p-value",
                       "adj Portm. stat", "adj Portm. p-value", 
                       "BG-LM stat", "BG LM-p-value")
  return(output)
}

####################################################################
#######    Long-run and short-run restriction matrix      ##########
####################################################################
#VARmodel = VAR6
responseMatrix = function(VARmodel) 
{
  Amats = Acoef(VARmodel)
  P     = VARmodel$p
  Ik    = diag(VARmodel$K)
  PsiMat  = matrix(0, VARmodel$K, VARmodel$K)
  
  for (i in 1:P) {
    PsiMat = PsiMat - Amats[[i]]
  }
  PsiMat  = Ik + PsiMat
  PsiMat  = solve(PsiMat)
  
  PsiMat[1,2:3] = 0
  PsiMat[2,3] = 0
  
  #n = 100 # slower
  #PsiMat = matrix(0,VARmodel$K,VARmodel$K)
  #for (i in 1:n) {
  #  PsiMat = PsiMat + Phi(x,nstep = i)[,,i]
  #}
  
  df             = summary(VARmodel$varresult[[1]])$df[2]
  SigmaU         = crossprod(resid(VARmodel))/df
  SigmaU         = crossprod(resid(VARmodel))/df
  LRIM           = solve(PsiMat) %*% SigmaU %*% solve(t(PsiMat))
  #LRIM           = PsiMat %*% SigmaU %*% t(PsiMat)
  colnames(LRIM) = colnames(x$y)
  rownames(LRIM) = colnames(lrim)
  SRIM           = PsiMat%*%LRIM           # still needs to be multiplied by restricted LRIM
  colnames(SRIM) = colnames(LRIM)
  rownames(SRIM) = colnames(LRIM)
  
  result         = list(A = Ik, B = SRIM, lrim = LRIM, Sigma.U = SigmaU * 100)
  return(result)
}