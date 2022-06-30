##################################
# Replication code for
# "Analyzing Attention to Scandal on Twitter"
#
# Last run on R version 4.1.1
#################################

# Contents:
# 0. load packages and data
# 1. Prepare time series
# 2. Run diagnostics on time series for modeling
# 3. Granger causality testing
# 4. Specify VARs
# 5. Collect results from cIRFs



###

# 0. Load packages and data

# Load packages

#install.packages('vars')
#install.packages('vctrs')
#install.packages('haven')
#install.packages('dplyr')
#install.packages('vars')
#install.packages('boot')
#install.packages('rio')

library(vars)
library(vctrs)
library(haven)
library(dplyr)
library(vars)
library(boot)
library(rio)

# Load custom functions

source("functions/granger.test.R")
source("functions/var.lag.specification.R")

# Load data, select variables of interest

db <- read.csv("timeseries.csv")
variables <- c("dem", "rep", "pubrep", "public", "media")
db <- db[, c("date", "topic", variables)]



###

# 1. Prepare scandals time series 

# Log daily attention

for (v in variables) {
  # - pulling the series-agenda for that group
  x <- db[,v]
  # - for some groups the last couple observations for each issues are NA,
  #     making these a 0 
  x[which(is.na(x))] <-0.01
  # - adding 1 percentage point to avoid 0s before the logit transformation
  #x <- x + 0.01
  # - applying the non-linear transformation
  logit_x <- log(x / (1-x))
  db[,v] <- logit_x
}



# Create scandal-specific datasets and time series

beng_data <- db[db$topic == 'benghazi',]
va_data <- db[db$topic == 'va',]
nsa_data <- db[db$topic == 'nsa',]
irs_data <- db[db$topic == 'irs',]

X_endog_beng <- beng_data[, which(grepl("rep|pubrep|media", colnames(db)))]
X_endog_va <- va_data[, which(grepl("rep|pubrep|media", colnames(db)))]
X_endog_nsa <- nsa_data[, which(grepl("rep|pubrep|media", colnames(db)))]
X_endog_irs <- irs_data[, which(grepl("rep|pubrep|media", colnames(db)))]



###

# 2. Run diagnostics on time series

# Phillips-Peron Ztau statistics

ur.pp(beng_data$rep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(beng_data$pubrep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(beng_data$media, type = "Z-tau", model = "constant", lags = "short")

ur.pp(va_data$rep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(va_data$pubrep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(va_data$media, type = "Z-tau", model = "constant", lags = "short")

ur.pp(nsa_data$rep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(nsa_data$pubrep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(nsa_data$media, type = "Z-tau", model = "constant", lags = "short")

ur.pp(irs_data$rep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(irs_data$pubrep, type = "Z-tau", model = "constant", lags = "short")
ur.pp(irs_data$media, type = "Z-tau", model = "constant", lags = "short")

# Determine lag lengths

VARselect(X_endog_beng, lag.max = 14, type = "const")
VARselect(X_endog_va, lag.max = 14, type = "const")
VARselect(X_endog_nsa, lag.max = 14, type = "const")
VARselect(X_endog_irs, lag.max = 14, type = "const")

# Specify VARs

setseed(2000)

var.beng <- VAR(y = X_endog_beng, p = 7, type="const")
var.va <- VAR(y = X_endog_va, p = 8, type="const") 
var.nsa <- VAR(y = X_endog_nsa, p = 7, type="const") 
var.irs <- VAR(y = X_endog_irs, p = 8, type="const")



###

# 3. Granger causality tests

granger.test(X_endog_beng,p=7)
granger.test(X_endog_va,p=8)
granger.test(X_endog_nsa,p=7)
granger.test(X_endog_irs,p=8)



###

# 4. Set up VARS with each actor ranked first by impulse status

rep.var.beng <- var.beng
rep.var.irs <- var.irs
rep.var.va <- var.va
rep.var.nsa <- var.nsa


pub.beng <- X_endog_beng[,c(2,1,3)]
pub.irs <- X_endog_irs[,c(2,1,3)]
pub.va <- X_endog_va[,c(2,1,3)]
pub.nsa <- X_endog_nsa[,c(2,1,3)]
pub.var.beng <- VAR(y = pub.beng, p = 7, type="const")
pub.var.irs <- VAR(y = pub.irs, p = 8, type="const")
pub.var.va <- VAR(y = pub.va, p = 7, type="const")
pub.var.nsa <- VAR(y = pub.nsa, p = 8, type="const")


media.beng <- X_endog_beng[,c(3,2,1)]
media.irs <- X_endog_irs[,c(3,2,1)]
media.va <- X_endog_va[,c(3,2,1)]
media.nsa <- X_endog_nsa[,c(3,2,1)]
media.var.beng <- VAR(y = media.beng, p = 7, type="const")
media.var.irs <- VAR(y = media.irs, p = 8, type="const")
media.var.va <- VAR(y = media.va, p = 7, type="const")
media.var.nsa <- VAR(y = media.nsa, p = 8, type="const")


###

# 5. Collect results from cIRFs

# prepare multiplier that transforms pulse of inv.logit(1) - or 0.73 - to 10 percentage points

multiplier <- (1/inv.logit(1))*10

# elites -> pubrep

beng1 <- irf(rep.var.beng, impulse="rep", response="pubrep", n.ahead = 15, cumulative = TRUE)
beng1 <- data.frame(cbind(beng1$irf$rep, beng1$Lower$rep, beng1$Upper$rep))
names(beng1) <- c("Value","lower","upper")
(inv.logit(beng1$Value)[15] - 0.5)*multiplier 
(inv.logit(beng1$upper)[15] - 0.5)*multiplier
(inv.logit(beng1$lower)[15] - 0.5)*multiplier


irs1 <- irf(rep.var.irs, impulse="rep", response="pubrep", n.ahead = 15, cumulative = TRUE)
irs1 <- data.frame(cbind(irs1$irf$rep, irs1$Lower$rep, irs1$Upper$rep))
names(irs1) <- c("Value","lower","upper")
inv.logit((irs1$Value)[15] - 0.5)*multiplier
inv.logit((irs1$upper)[15] - 0.5)*multiplier
inv.logit((irs1$lower)[15] - 0.5)*multiplier

va1 <- irf(rep.var.va, impulse="rep", response="pubrep", n.ahead = 15, cumulative = TRUE)
va1 <- data.frame(cbind(va1$irf$rep, va1$Lower$rep, va1$Upper$rep))
names(va1) <- c("Value","lower","upper")
(inv.logit(va1$Value)[15] - 0.5)*multiplier
(inv.logit(va1$upper)[15] - 0.5)*multiplier
(inv.logit(va1$lower)[15] - 0.5)*multiplier

nsa1 <- irf(rep.var.nsa, impulse="rep", response="pubrep", n.ahead = 15, cumulative = TRUE)
nsa1 <- data.frame(cbind(nsa1$irf$rep, nsa1$Lower$rep, nsa1$Upper$rep))
names(nsa1) <- c("Value","lower","upper")
inv.logit((nsa1$Value)[15] - 0.5)*multiplier
inv.logit((nsa1$upper)[15] - 0.5)*multiplier
inv.logit((nsa1$lower)[15] - 0.5)*multiplier

# pubrep -> elites

beng2.1 <- irf(pub.var.beng, impulse="pubrep", response="rep", n.ahead = 15, cumulative = TRUE)
beng2.1 <- data.frame(cbind(beng2.1$irf$pubrep, beng2.1$Lower$pubrep, beng2.1$Upper$pubrep))
names(beng2.1) <- c("Value","lower","upper")
(inv.logit(beng2.1$Value) - 0.5)[15]*multiplier
(inv.logit(beng2.1$upper) - 0.5)[15]*multiplier
(inv.logit(beng2.1$lower) - 0.5)[15]*multiplier

irs2.1 <- irf(pub.var.irs, impulse="pubrep", response="rep", n.ahead = 15, cumulative = TRUE)
irs2.1 <- data.frame(cbind(irs2.1$irf$pubrep, irs2.1$Lower$pubrep, irs2.1$Upper$pubrep))
names(irs2.1) <- c("Value","lower","upper")
(inv.logit(irs2.1$Value) - 0.5)[15]*multiplier
(inv.logit(irs2.1$upper) - 0.5)[15]*multiplier
(inv.logit(irs2.1$lower) - 0.5)[15]*multiplier

va2.1 <- irf(pub.var.va, impulse="pubrep", response="rep", n.ahead = 15, cumulative = TRUE)
va2.1 <- data.frame(cbind(va2.1$irf$pubrep, va2.1$Lower$pubrep, va2.1$Upper$pubrep))
names(va2.1) <- c("Value","lower","upper")
(inv.logit(va2.1$Value) - 0.5)[15]*multiplier
(inv.logit(va2.1$upper) - 0.5)[15]*multiplier
(inv.logit(va2.1$lower) - 0.5)[15]*multiplier

nsa2.1 <- irf(pub.var.nsa, impulse="pubrep", response="rep", n.ahead = 15, cumulative = TRUE)
nsa2.1 <- data.frame(cbind(nsa2.1$irf$pubrep, nsa2.1$Lower$pubrep, nsa2.1$Upper$pubrep))
names(nsa2.1) <- c("Value","lower","upper")
(inv.logit(nsa2.1$Value) - 0.5)[15]*multiplier
(inv.logit(nsa2.1$upper) - 0.5)[15]*multiplier
(inv.logit(nsa2.1$lower) - 0.5)[15]*multiplier


# pubrep -> media

beng2.2 <- irf(pub.var.beng, impulse="pubrep", response="media", n.ahead = 15, cumulative = TRUE)
beng2.2 <- data.frame(cbind(beng2.2$irf$pubrep, beng2.2$Lower$pubrep, beng2.2$Upper$pubrep))
names(beng2.2) <- c("Value","lower","upper")
(inv.logit(beng2.2$Value) - 0.5)[15]*multiplier
(inv.logit(beng2.2$upper) - 0.5)[15]*multiplier
(inv.logit(beng2.2$lower) - 0.5)[15]*multiplier

irs2.2 <- irf(pub.var.irs, impulse="pubrep", response="media", n.ahead = 15, cumulative = TRUE)
irs2.2 <- data.frame(cbind(irs2.2$irf$pubrep, irs2.2$Lower$pubrep, irs2.2$Upper$pubrep))
names(irs2.2) <- c("Value","lower","upper")
(inv.logit(irs2.2$Value) - 0.5)[15]*multiplier
(inv.logit(irs2.2$upper) - 0.5)[15]*multiplier
(inv.logit(irs2.2$lower) - 0.5)[15]*multiplier

va2.2 <- irf(pub.var.va, impulse="pubrep", response="media", n.ahead = 15, cumulative = TRUE)
va2.2 <- data.frame(cbind(va2.2$irf$pubrep, va2.2$Lower$pubrep, va2.2$Upper$pubrep))
names(va2.2) <- c("Value","lower","upper")
(inv.logit(va2.2$Value) - 0.5)[15]*multiplier
(inv.logit(va2.2$upper) - 0.5)[15]*multiplier
(inv.logit(va2.2$lower) - 0.5)[15]*multiplier

nsa2.2 <- irf(pub.var.nsa, impulse="pubrep", response="media", n.ahead = 15, cumulative = TRUE)
nsa2.2 <- data.frame(cbind(nsa2.2$irf$pubrep, nsa2.2$Lower$pubrep, nsa2.2$Upper$pubrep))
names(nsa2.2) <- c("Value","lower","upper")
(inv.logit(nsa2.2$Value) - 0.5)[15]*multiplier
(inv.logit(nsa2.2$upper) - 0.5)[15]*multiplier
(inv.logit(nsa2.2$lower) - 0.5)[15]*multiplier

# elites -> media

beng3.2 <- irf(rep.var.beng, impulse="rep", response="media", n.ahead = 15, cumulative = TRUE)
beng3.2 <- data.frame(cbind(beng3.2$irf$rep, beng3.2$Lower$rep, beng3.2$Upper$rep))
names(beng3.2) <- c("Value","lower","upper")
(inv.logit(beng3.2$Value) - 0.5)[15]*multiplier
(inv.logit(beng3.2$upper) - 0.5)[15]*multiplier
(inv.logit(beng3.2$lower) - 0.5)[15]*multiplier

irs3.2 <- irf(rep.var.irs, impulse="rep", response="media", n.ahead = 15, cumulative = TRUE)
irs3.2 <- data.frame(cbind(irs3.2$irf$rep, irs3.2$Lower$rep, irs3.2$Upper$rep))
names(irs3.2) <- c("Value","lower","upper")
(inv.logit(irs3.2$Value) - 0.5)[15]*multiplier
(inv.logit(irs3.2$upper) - 0.5)[15]*multiplier
(inv.logit(irs3.2$lower) - 0.5)[15]*multiplier

va3.2 <- irf(rep.var.va, impulse="rep", response="media", n.ahead = 15, cumulative = TRUE)
va3.2 <- data.frame(cbind(va3.2$irf$rep, va3.2$Lower$rep, va3.2$Upper$rep))
names(va3.2) <- c("Value","lower","upper")
(inv.logit(va3.2$Value) - 0.5)[15]*multiplier
(inv.logit(va3.2$upper) - 0.5)[15]*multiplier
(inv.logit(va3.2$lower) - 0.5)[15]*multiplier

nsa3.2 <- irf(rep.var.nsa, impulse="rep", response="media", n.ahead = 15, cumulative = TRUE)
nsa3.2 <- data.frame(cbind(nsa3.2$irf$rep, nsa3.2$Lower$rep, nsa3.2$Upper$rep))
names(nsa3.2) <- c("Value","lower","upper")
(inv.logit(nsa3.2$Value) - 0.5)[15]*multiplier
(inv.logit(nsa3.2$upper) - 0.5)[15]*multiplier
(inv.logit(nsa3.2$lower) - 0.5)[15]*multiplier

# media -> pubrep

beng4 <- irf(media.var.beng, impulse="media", response="pubrep", n.ahead = 15, cumulative = TRUE)
beng4 <- data.frame(cbind(beng4$irf$media, beng4$Lower$media, beng4$Upper$media))
names(beng4) <- c("Value","lower","upper")
(inv.logit(beng4$Value) - 0.5)[15]*multiplier
(inv.logit(beng4$upper) - 0.5)[15]*multiplier
(inv.logit(beng4$lower) - 0.5)[15]*multiplier

irs4 <- irf(media.var.irs, impulse="media", response="pubrep", n.ahead = 15, cumulative = TRUE)
irs4 <- data.frame(cbind(irs4$irf$media, irs4$Lower$media, irs4$Upper$media))
names(irs4) <- c("Value","lower","upper")
(inv.logit(irs4$Value) - 0.5)[15]*multiplier
(inv.logit(irs4$upper) - 0.5)[15]*multiplier
(inv.logit(irs4$lower) - 0.5)[15]*multiplier

va4 <- irf(media.var.va, impulse="media", response="pubrep", n.ahead = 15, cumulative = TRUE)
va4 <- data.frame(cbind(va4$irf$media, va4$Lower$media, va4$Upper$media))
names(va4) <- c("Value","lower","upper")
(inv.logit(va4$Value) - 0.5)[15]*multiplier
(inv.logit(va4$upper) - 0.5)[15]*multiplier
(inv.logit(va4$lower) - 0.5)[15]*multiplier

nsa4 <- irf(media.var.nsa, impulse="media", response="pubrep", n.ahead = 15, cumulative = TRUE)
nsa4 <- data.frame(cbind(nsa4$irf$media, nsa4$Lower$media, nsa4$Upper$media))
names(nsa4) <- c("Value","lower","upper")
(inv.logit(nsa4$Value) - 0.5)[15]*multiplier
(inv.logit(nsa4$upper) - 0.5)[15]*multiplier
(inv.logit(nsa4$lower) - 0.5)[15]*multiplier

# media -> elites


beng5 <- irf(media.var.beng, impulse="media", response="rep", n.ahead = 15, cumulative = TRUE)
beng5 <- data.frame(cbind(beng5$irf$media, beng5$Lower$media, beng5$Upper$media))
names(beng5) <- c("Value","lower","upper")
(inv.logit(beng5$Value) - 0.5)[15]*multiplier
(inv.logit(beng5$upper) - 0.5)[15]*multiplier
(inv.logit(beng5$lower) - 0.5)[15]*multiplier

irs5 <- irf(media.var.irs, impulse="media", response="rep", n.ahead = 15, cumulative = TRUE)
irs5 <- data.frame(cbind(irs5$irf$media, irs5$Lower$media, irs5$Upper$media))
names(irs5) <- c("Value","lower","upper")
(inv.logit(irs5$Value) - 0.5)[15]*multiplier
(inv.logit(irs5$upper) - 0.5)[15]*multiplier
(inv.logit(irs5$lower) - 0.5)[15]*multiplier

va5 <- irf(media.var.va, impulse="media", response="rep", n.ahead = 15, cumulative = TRUE)
va5 <- data.frame(cbind(va5$irf$media, va5$Lower$media, va5$Upper$media))
names(va5) <- c("Value","lower","upper")
(inv.logit(va5$Value) - 0.5)[15]*multiplier
(inv.logit(va5$upper) - 0.5)[15]*multiplier
(inv.logit(va5$lower) - 0.5)[15]*multiplier

nsa5 <- irf(media.var.nsa, impulse="media", response="rep", n.ahead = 15, cumulative = TRUE)
nsa5 <- data.frame(cbind(nsa5$irf$media, nsa5$Lower$media, nsa5$Upper$media))
names(nsa5) <- c("Value","lower","upper")
(inv.logit(nsa5$Value) - 0.5)[15]*multiplier
(inv.logit(nsa5$upper) - 0.5)[15]*multiplier
(inv.logit(nsa5$lower) - 0.5)[15]*multiplier
