library('rmgarch')
library('zoo')
library("forecast")
library('rugarch')
library('lmtest')
library('lattice')
library('MASS')
library('fitdistrplus')
library('metRology')
library('ggplot2')
library('reshape2')
library('parallel')
library('reticulate')
library('jsonlite')
library(stargazer)
library(urca)
library('car')
        
setwd("~/Desktop/Dissertation")



# Urca test interpretation ------- 
interp_urdf <- function(urdf, level="5pct") {
  if(class(urdf) != "ur.df") stop('parameter is not of class ur.df from urca package')
  if(!(level %in% c("1pct", "5pct", "10pct") ) ) stop('parameter level is not one of 1pct, 5pct, or 10pct')
  
  cat("========================================================================\n")
  cat( paste("At the", level, "level:\n") )
  if(urdf@model == "none") {
    cat("The model is of type none\n")
    tau1_crit = urdf@cval["tau1",level]
    tau1_teststat = urdf@teststat["statistic","tau1"]
    tau1_teststat_wi_crit = tau1_teststat > tau1_crit
    if(tau1_teststat_wi_crit) {
      cat("tau1: The null hypothesis is not rejected, unit root is present\n")
    } else {
      cat("tau1: The null hypothesis is rejected, unit root is not present\n")
    }
  } else if(urdf@model == "drift") {
    cat("The model is of type drift\n")
    tau2_crit = urdf@cval["tau2",level]
    tau2_teststat = urdf@teststat["statistic","tau2"]
    tau2_teststat_wi_crit = tau2_teststat > tau2_crit
    phi1_crit = urdf@cval["phi1",level]
    phi1_teststat = urdf@teststat["statistic","phi1"]
    phi1_teststat_wi_crit = phi1_teststat < phi1_crit
    if(tau2_teststat_wi_crit) {
      # Unit root present branch
      cat("tau2: The first null hypothesis is not rejected, unit root is present\n")
      if(phi1_teststat_wi_crit) {
        cat("phi1: The second null hypothesis is not rejected, unit root is present\n")
        cat("      and there is no drift.\n")
      } else {
        cat("phi1: The second null hypothesis is rejected, unit root is present\n")
        cat("      and there is drift.\n")
      }
    } else {
      # Unit root not present branch
      cat("tau2: The first null hypothesis is rejected, unit root is not present\n")
      if(phi1_teststat_wi_crit) {
        cat("phi1: The second null hypothesis is not rejected, unit root is present\n")
        cat("      and there is no drift.\n")
        warning("This is inconsistent with the first null hypothesis.")
      } else {
        cat("phi1: The second null hypothesis is rejected, unit root is not present\n")
        cat("      and there is drift.\n")
      }
    }
  } else if(urdf@model == "trend") {
    cat("The model is of type trend\n")
    tau3_crit = urdf@cval["tau3",level]
    tau3_teststat = urdf@teststat["statistic","tau3"]
    tau3_teststat_wi_crit = tau3_teststat > tau3_crit
    phi2_crit = urdf@cval["phi2",level]
    phi2_teststat = urdf@teststat["statistic","phi2"]
    phi2_teststat_wi_crit = phi2_teststat < phi2_crit
    phi3_crit = urdf@cval["phi3",level]
    phi3_teststat = urdf@teststat["statistic","phi3"]
    phi3_teststat_wi_crit = phi3_teststat < phi3_crit
    if(tau3_teststat_wi_crit) {
      # First null hypothesis is not rejected, Unit root present branch
      cat("tau3: The first null hypothesis is not rejected, unit root is present\n")
      if(phi3_teststat_wi_crit) {
        # Second null hypothesis is not rejected
        cat("phi3: The second null hypothesis is not rejected, unit root is present\n")
        cat("      and there is no trend\n")
        if(phi2_teststat_wi_crit) {
          # Third null hypothesis is not rejected
          # a0-drift = gamma = a2-trend = 0
          cat("phi2: The third null hypothesis is not rejected, unit root is present\n")
          cat("      there is no trend, and there is no drift\n")
        } else {
          # Third null hypothesis is rejected
          cat("phi2: The third null hypothesis is rejected, unit root is present\n")
          cat("      there is no trend, and there is drift\n")
        }
      }
      else {
        # Second null hypothesis is rejected
        cat("phi3: The second null hypothesis is rejected, unit root is present\n")
        cat("      and there is trend\n")
        if(phi2_teststat_wi_crit) {
          # Third null hypothesis is not rejected
          # a0-drift = gamma = a2-trend = 0
          cat("phi2: The third null hypothesis is not rejected, unit root is present\n")
          cat("      there is no trend, and there is no drift\n")
          warning("This is inconsistent with the second null hypothesis.")
        } else {
          # Third null hypothesis is rejected
          cat("phi2: The third null hypothesis is rejected, unit root is present\n")
          cat("      there is trend, and there may or may not be drift\n")
          warning("Presence of drift is inconclusive.")
        }
      }
    } else {
      # First null hypothesis is rejected, Unit root not present branch
      cat("tau3: The first null hypothesis is rejected, unit root is not present\n")
      if(phi3_teststat_wi_crit) {
        cat("phi3: The second null hypothesis is not rejected, unit root is present\n")
        cat("      and there is no trend\n")
        warning("This is inconsistent with the first null hypothesis.")
        if(phi2_teststat_wi_crit) {
          # Third null hypothesis is not rejected
          # a0-drift = gamma = a2-trend = 0
          cat("phi2: The third null hypothesis is not rejected, unit root is present\n")
          cat("      there is no trend, and there is no drift\n")
          warning("This is inconsistent with the first null hypothesis.")
        } else {
          # Third null hypothesis is rejected
          cat("phi2: The third null hypothesis is rejected, unit root is not present\n")
          cat("      there is no trend, and there is drift\n")
        }
      } else {
        cat("phi3: The second null hypothesis is rejected, unit root is not present\n")
        cat("      and there may or may not be trend\n")
        warning("Presence of trend is inconclusive.")
        if(phi2_teststat_wi_crit) {
          # Third null hypothesis is not rejected
          # a0-drift = gamma = a2-trend = 0
          cat("phi2: The third null hypothesis is not rejected, unit root is present\n")
          cat("      there is no trend, and there is no drift\n")
          warning("This is inconsistent with the first and second null hypothesis.")
        } else {
          # Third null hypothesis is rejected
          cat("phi2: The third null hypothesis is rejected, unit root is not present\n")
          cat("      there may or may not be trend, and there may or may not be drift\n")
          warning("Presence of trend and drift is inconclusive.")
        }
      }
    }
  } else warning('urdf model type is not one of none, drift, or trend')
  cat("========================================================================\n")
}






# Reading and transforming data -------------------------------------------
df <- data.frame(read.csv('data_stocks.csv', header=T))


# Converting data types of columns df.
numeric_columns <- seq(2, dim(df)[2])
df[,numeric_columns] <- lapply(df[,numeric_columns], function(x) as.numeric(x))
df[,-numeric_columns] <- as.Date(df$Date, format="%Y-%m-%d")

# Train, test datasets.
df_train <- df[1:(0.80*length(df$Date)),]
indices_train <- as.numeric(row.names(df_train))
df_test <- df[-indices_train, ]

# Creating data.frame of log changes.
df_linear_returns <- apply(df[,numeric_columns], 2, function(x)  ((as.numeric(x[-1]) - as.numeric(x[-length(x)]))/as.numeric(x[-length(x)])))
df_linear_returns <- data.frame(subset(df, select='Date', row.names(df) > 1), df_linear_returns)

df_linear_returns_train <- df_linear_returns[indices_train, ]
df_linear_returns_test <- df_linear_returns[-(indices_train), ]

# Plotting time-series AMZN.
zoo_amzn <- zoo(df_linear_returns$AMZN, df_linear_returns$Date)
plot(zoo_amzn, xlab='Time', ylab='Linear returns of AMZN')

# ARIMA model -------------------------------------------------------------
# ACF and PACF:

par(mfrow=c(1,1))
Acf(df_linear_returns_train$AMZN, main='ACF of AMZN')
pacf(df_linear_returns_train$AMZN, main='PACF of AMZN')


# ARIMA model.
model1 <- Arima(df_linear_returns_train$AMZN, include.mean = T, include.drift = T)
summary(model1)
tsdisplay(model1$residuals)

# ADF-test:
adf_test1 <- (ur.df(df_linear_returns$FANG, type = "trend"))
interp_urdf(adf_test1)
summary(adf_test1)
# ADF test:
adf_test2 <- (ur.df(diff(df_linear_returns$AMZN), type = "trend"))
interp_urdf(adf_test2)


# Auto-ARIMA model.
model_auto <- auto.arima(df_linear_returns_train$AMZN, stepwise = FALSE, parallel = T, test = c("adf"),
                                ic = c("aicc"))

model_auto
# Looking for the heteroscedasticity.
par(mfrow=c(1,2))
Acf(model1$residuals^2, lag=50, main='ACF of square of residuals from ARMA model')
Acf(abs(model1$residuals), lag=50, main='ACF of square of residuals from ARMA model')
par(mfrow=c(1,1))



# Fitting distribution ----------------------------------------------------
# Fitting data to the t distribution.
fit.distr <- fitdist(as.numeric(model1$residuals), "t.scaled", method = 'mle', start = list(df = 2, mean = 1, sd = 1))
plot(fit.distr, breaks = 50)
qqPlot(as.numeric(model1$residuals), distribution = 't.scaled', df = fit.distr$estimate['df'], 
       mean = fit.distr$estimate['mean'], sd = fit.distr$estimate['sd'])
# Theoretical density plot and histogram.
res <- as.numeric(model1$residuals)
theoretical_x <- rt.scaled(10000, df = fit.distr$estimate['df'], 
                          mean = fit.distr$estimate['mean'], sd = fit.distr$estimate['sd'])

h <- hist(as.numeric(res), freq=FALSE, breaks = 50, prob = TRUE)
lines(density(theoretical_x))


# ARMA + GARCH FANG------------------------------------------------------------

cl = makePSOCKcluster(10)

ctrl<- list(RHO = 1, DELTA = 1e-8, MAJIT = 100, MINIT = 650, TOL = 1e-6)

garch.spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                               garchOrder = c(1,1)), 
                         mean.model=list(armaOrder = c(0,0), 
                                         include.mean = TRUE), 
                         distribution.model = "std") 

garch.fit <- ugarchfit(data = zoo(df_linear_returns_train$AMZN, df_linear_returns_train$Date), 
                       spec = garch.spec, 
                       solver = "solnp", 
                       solver.control = ctrl, cluster = cl)  

garch.fit
## QQ-plot of standardized residuals
plot(garch.fit, which = 9)
## ACF plot of standardized residuals(for unconditional heteroscedasticity).
plot(garch.fit, which = 10)
## ACF plot of squared standardized residuals(for conditional heteroscedasticity)
plot(garch.fit, which = 11)

# Predictions for several 
forecast_length <- as.numeric(difftime(as.POSIXct(max(df_test$Date)), as.POSIXct(min(df_test$Date)), units = 'days'))
forc.res <- ugarchforecast(garch.fit, n.ahead = forecast_length)
## take it as a data frame
forc<- as.data.frame(forc.res@forecast)
names(forc)<-c("n.ahead","N","n.start","n.roll","sigma","series")

# Forecasting ARMA + GARCH model-------------------------------------------------------------
linear_transform <- function(forecast_length, original_data, linear_data, name_of_column){
  transformed_data = tail(original_data[, name_of_column], 1)
  for(i in seq(forecast_length)){
    transformed_data = append(transformed_data, (linear_data[i] * transformed_data[i] + transformed_data[i]))
  }
  transformed_data <- tail(transformed_data, -1)
  return(transformed_data)
}

price_forecast <- linear_transform(forecast_length, df_train, forc$series, 'AMZN')
# Accuracy of forecast:
accuracy(price_forecast, df_test$AMZN)

# Plotting of forecast.
par(mfrow=c(1,1))
eqix.index <- df_train$Date
eqix.values <- df_train$AMZN
eqix <- zoo(eqix.values, eqix.index)
eqix_future.index <- df_test$Date
eqix_future.values <- df_test$AMZN
eqix_future <- zoo(eqix_future.values, eqix_future.index)
eqix_forecast <- zooreg(price_forecast, start  = as.Date(tail(df_train$Date, 1)))
plot.zoo(cbind(eqix, eqix_future), col = c('black', 'blue'), plot.type = 'single', main='1 year ahead forecast for AMZN',
         )
lines(eqix_forecast, col = 'red')


# Generating scenarios based on simulation.
# Since we make simulation between tail(df_train$Date, 1) and (head(df_test$Date, 1) additional simulatio.
sim_10_log <- ugarchsim(garch.fit, n.sim = forecast_length + 1, m.sim = 10)
sim_10 <- apply(sim_10_log@simulation$seriesSim, MARGIN = 2, linear_transform, forecast_length = forecast_length + 1, 
       original_data = df_train, name_of_column = "AMZN")

# Transforming into dataframe.
seq_dates <- seq(min(df_test$Date), max(df_test$Date), by = 1)
data_sim_10 <- data.frame(simulation = sim_10, Date = seq_dates)
head(data_sim_10)

# Plotting path of scenarios.
plot.zoo(eqix, xlim = c(min(df$Date), max(df$Date)), ylim = (c(min(df$AMZN), 800)))
for (i in seq(1, 10)){
  lines(zooreg(sim_10[, i], start = as.Date(tail(df_train$Date, 1))), col = i + 2)
}


meltdf <- melt(data_sim_10, id="Date")
ggplot() + 
  geom_line(data = data.frame(AMZN = eqix.values, Date = eqix.index), mapping = aes(x=Date, y=AMZN)) +
  geom_line(data = meltdf, aes(x=Date, y=value, colour=variable, group=variable))



# DCC-GARCH ---------------------------------------------------------------
garch.spec2 <- ugarchspec(variance.model = list(model = "sGARCH", 
                                               garchOrder = c(1,1)), 
                         mean.model=list(armaOrder = c(0,0), 
                                         include.mean = TRUE), 
                         distribution.model = "std",) 

garch.spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                               garchOrder = c(1,1)), 
                         mean.model=list(armaOrder = c(0,0), 
                                         include.mean = TRUE), 
                         distribution.model = "std") 

uspec.n2 = multispec(replicate(6, garch.spec2))
uspec.n = multispec(replicate(6, garch.spec))

spec.dccn2 = dccspec(uspec.n2, dccOrder = c(1, 1), distribution = 'mvt', model = 'DCC')
spec.dccn = dccspec(uspec.n, dccOrder = c(1, 1), distribution = 'mvt', model = 'DCC')

# Transforming df_linear_returns into zoo object. 
zoo_linear_returns <- zoo(x = df_linear_returns_train[,2:dim(df_linear_returns_train)[2]], order.by = df_linear_returns_train$Date)
zoo_linear <- zoo(df_linear_returns[,-1], order.by = df_linear_returns$Date)

dcc_model <- dccfit(spec.dccn, data = zoo_linear_returns, solver = 'solnp', cluster = cl)

plot(dcc_model, which = 5, series=seq(1, 6))

# Forecasting.

dcc_forecast <- dccforecast(dcc_model, n.ahead = forecast_length)
dcc_forecast <- dcc_forecast@mforecast$mu

linear_multi_transform <- function(forecast_length, original_data, linear_data){
  transformed_data = original_data[dim(original_data)[1], -1]
  row.names(transformed_data) <- NULL
  for(i in seq(forecast_length)){
    transformed_data = rbind(transformed_data, (linear_data[i, ,1] * transformed_data[i, ] + transformed_data[i, ]))
    row.names(transformed_data) <- NULL
  }
  transformed_data <- tail(transformed_data, -1)
  return(transformed_data)
}

price_forecast <- linear_multi_transform(forecast_length, df_train, dcc_forecast)

# Plotting the forecast:

train_set <- zoo(df_train[, -1], order.by = df_train$Date)
test_set <- zoo(df_test[, -1], order.by = df_test$Date)
forecast_set <- zooreg(price_forecast, start  = as.Date(tail(df_train$Date, 1)))
par(mfrow=c(2,3))
for(i in seq(dim(forecast_set)[2])){
  plot(cbind(train_set[, i], test_set[, i]), col = c('black', 'blue'), plot.type = 'single', 
       main = colnames(df_train[, -1])[i], xlab = 'Date', ylab = 'Price')
  lines(forecast_set[, i], col = 'red')
  }


par(mfrow=c(1,1))

# Student Copula - GARCH ---------------------------------------------------------------


spec.copula = cgarchspec(uspec.n2, dccOrder = c(1, 1), distribution.model = list(copula = 'mvt'))
copula_model <- cgarchfit(spec.copula, data = zoo_linear_returns, cluster = cl, fit.control = list(eval.se = FALSE))
cgarchsim(copula_model, n.sim = 1490, m.sim = 10, n.start = 0, startMethod = 'sample', rseed = 1)


# Roll forecast.
dcc_roll <- dccroll(spec.dccn, data = zoo_linear, refit.every = 30, n.start = (dim(df_linear_returns_train)[1] - 1), n.ahead = 1,
        refit.window = 'moving', window.size = 300, cluster = cl, save.fit = FALSE, clusterOnAssets = TRUE)

dcc_roll2 <- dccroll(spec.dccn2, data = zoo_linear, refit.every = 30, n.start = (dim(df_linear_returns_train)[1] - 1), n.ahead = 1,
        refit.window = 'moving', window.size = 300, cluster = cl, save.fit = TRUE, save.wdir = 'results_dcc', clusterOnAssets = TRUE)

dcc_roll_forecast <- data.frame(coredata(fitted(dcc_roll)))
dcc_roll_sigma <- coredata(rcov(dcc_roll))
dcc_sigma <- rep(list(), dim(dcc_roll_sigma)[3])
for(i in seq(dim(dcc_roll_sigma)[3])){
  dcc_sigma[[i]] <- data.frame(dcc_roll_sigma[, , i])
}

# Writing to JSON data from dcc roll.
# write_json(dcc_roll_forecast, "set_sim_mean.json", digits = 7)
# write_json(dcc_sigma, 'set_sim_covar.json', digits = 7, dataframe = 'column')
# write_json(df_linear_returns_test, 'df_test.json', digits = 7)


# Saving refit models.
list_files <- list.files(path = "results_dcc")
set_dcc_models <- c()
for(i in seq(length(list_files))){
  set_dcc_models <- c(set_dcc_models, readRDS(paste('results_dcc/', list_files[i], sep = '')))
}


# Plotting DCC roll forecast ----------------------------------------------

train_set <- zoo(df_train[, -1], order.by = df_train$Date)
test_set <- zoo(df_test[, -1], order.by = df_test$Date)

linear_multi_transform <- function(forecast_length, original_data, linear_data, realized_data){
  data <- original_data[dim(original_data)[1], -1]
  forecast_data <- data
  data <- rbind(data, realized_data[,-1])
  row.names(forecast_data) <- NULL
  row.names(data) <- NULL
  for(i in seq(forecast_length)){
    forecast_data <- rbind(forecast_data, (linear_data[i, ] * data[i, ] + data[i, ]))
    row.names(forecast_data) <- NULL
  }
  forecast_data <- tail(forecast_data, -1)
  return(forecast_data)
}

price_forecast <- linear_multi_transform(dim(dcc_roll_forecast)[1], df_train, data.frame(dcc_roll_forecast), df_test[])

forecast_set <- zoo(price_forecast, order.by = df_test$Date)
par(mfrow=c(2,3))
for(i in seq(dim(forecast_set)[2])){
  plot(cbind(train_set[, i], test_set[, i]), col = c('black', 'blue'), plot.type = 'single', 
       main = colnames(df_train[, -1])[i], xlab = 'Date', ylab = 'Price')
  lines(forecast_set[, i], col = 'red')
  }
title(sprintf("refit.every = %s, refit.window = %s, window.size = %s", dcc_roll@model[["refit.every"]], 
              dcc_roll@model[["refit.window"]], dcc_roll@model[["window.size"]]), outer = TRUE, line = -1)

# Plotting difference between one ahead forecast and realized price. 
# Vertical difference
plot(zoo(df_test[, -1] - price_forecast, order.by = df_test$Date), xlab='Date', main = "Difference between one ahead forecast and realized prices")

# On the same plot.
par(mfrow=c(2,3))
for (i in seq(1, 6)){
  plot(test_set[,i], col = 'blue')
  lines(forecast_set[,i], col = 'red')
}
par(mfow=c(1,1))

# One-ahead simulating rolling DCC-GARCH 
zoo_linear_for_roll <- head(c(tail(zoo_linear_returns, 1), zoo_linear[-seq(1, dim(zoo_linear_returns)[1])]), -1)
sim_number <- 200
for(i in seq(length(set_dcc_models))){
  n.start <- set_dcc_models[[i]]@model[["modeldata"]][["T"]]
  sim_horizon <- set_dcc_models[[i]]@model[["modeldata"]][["n.start"]]
  dcc_sim <- dccsim(set_dcc_models[[i]], n.sim = sim_horizon, m.sim = sim_number, n.start = n.start, startMethod = "sample", rseed = 123,
                    cluster = cl)}
dcc_sim@msim$simX[[28]]
# Roll forecast for n.ahead > 1 -------------------------------------------

sim_number <- 200
set_sim_mean <- rep(list(data.frame()), length(df_test$Date))
set_sim_covar <- rep(list(list()), length(df_test$Date))
len_test <- length(df_test$Date)
set_sim_covar <- list()
m <- 1

for(i in seq(1, length(set_dcc_models))){
  n.start <- set_dcc_models[[i]]@model[["modeldata"]][["T"]]
  sim_horizon <- set_dcc_models[[i]]@model[["modeldata"]][["n.start"]]
  j <- sim_horizon
  
  while(j >= 1){
    presigma <- sigma(set_dcc_models[[i]])[n.start  - j, ]
    preresiduals <- residuals(set_dcc_models[[i]])[n.start - j, ]
    prereturns <- matrix(set_dcc_models[[i]]@model$modeldata$data[n.start + sim_horizon - j, ], nrow = 1, ncol = 6, byrow = TRUE)
    preQ <- set_dcc_models[[i]]@mfit[["Q"]][[n.start - j]]
    Qbar <- set_dcc_models[[i]]@mfit$Qbar
    preZ <- set_dcc_models[[i]]@mfit$stdresid[n.start - j, ]
    
    examp_sim <- dccsim(set_dcc_models[[i]], n.sim = 1, m.sim = sim_number, startMethod = "sample", n.start = 1,
                        rseed = 1, cluster = cl, presigma = presigma, preresiduals = preresiduals, 
                        prereturns = prereturns, preQ = preQ, Qbar = Qbar, preZ = preZ)
    
    j <- j - 1
    set_sim_covar[[m]] <- rep(list(list()), sim_horizon)
    for(l in seq(sim_number)){
      set_sim_mean[[m]] <- na.omit(rbind(set_sim_mean[[m]], examp_sim@msim$simX[[l]]))
      set_sim_covar[[m]][[l]] <- examp_sim@msim$simH[[l]][, ,1]
    }
    row.names(set_sim_mean[[m]]) <- NULL
    colnames(set_sim_mean[[m]]) <- colnames(df[,-1])
    m <- m + 1
    }
}



# Transefering data to python ---------------------------------------------

# write_json(set_sim_mean, "set_montecarlo_mean.json", digits = 7)
# write_json(set_sim_covar[-length(set_sim_covar)], 'set_montecarlo_covar.json', digits = 7)
# write_json(df_linear_returns_test, 'df_test.json', digits = 7)



# Experiments -------
data(dji30retw)
Dat = dji30retw[, 1:3, drop = FALSE]
cnames = colnames(Dat)
uspec = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = T), 
                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
                   distribution.model = "std")
spec1 = dccspec(uspec = multispec( replicate(3, uspec) ), dccOrder = c(1,1), 
                distribution = "mvt")
fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se=FALSE))

presigma = tail( sigma(fit1 ), 1 )
preresiduals = tail( residuals(fit1), 1 )
prereturns = tail( as.matrix(Dat), 1 )
sim1 = dccsim(fitORspec = fit1, n.sim = 1000, n.start = 100, m.sim = 2, startMethod = "unconditional", 
              presigma = presigma, preresiduals = preresiduals, prereturns = prereturns, 
              preQ = last(rcor(fit1, type = "Q"))[,,1], Qbar = fit1@mfit$Qbar, 
              preZ = tail(fit1@mfit$stdresid, 1),
              rseed = c(100, 200), mexsimdata = NULL, vexsimdata = NULL)


ugarchspec(variance.model = list(model = "sGARCH", 
                                 garchOrder = c(1,1)), 
           mean.model=list(armaOrder = c(1,1), 
                           include.mean = TRUE, arfima = T), 
           distribution.model = "std", fixed.pars = list(arfima = 1)) 
