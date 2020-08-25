library(hdm)
library(restorepoint)
source("lasso_tools.R")
source("orthoML.R")


data(GrowthData)
dim(GrowthData)
## [1] 90 63
y = GrowthData[, 1, drop = T]
d = GrowthData[, 3, drop = T]
X = as.matrix(GrowthData)[, -c(1, 2, 3)]
dX = as.matrix(GrowthData)[, -c(1, 2)]

varnames = colnames(GrowthData)
dvar = varnames[3]
xnames = varnames[-c(1, 2, 3)] # names of X variables
dandxnames = varnames[-c(1, 2)] # names of D and X variables


# create formulas by pasting names (this saves typing times)
fmla = as.formula(paste("Outcome ~ ", paste(dandxnames, collapse = "+")))
ls.effect = lm(fmla, data = GrowthData)
coef(ls.effect)[2]

# Double selection
lasso.effect = rlassoEffect(x = X, y = y, d = d, method = "double selection")
summary(lasso.effect) # 0.05001

# Partialling out
lasso.effect = rlassoEffect(x = X, y = y, d = d, method = "partialling out")
summary(lasso.effect) # 0.04981

# Normal lasso using gamlr
lasso.gamlr = gamlr(x=dX,y=y,free=1)
post_lasso_coef(lasso.gamlr, dX,y)[1] # -0.05023

# Taddy's orthoML
set.seed(1)
orthoML( x=X, d=d, y=y, nfold=5)$coef # -0.0657


# As a test custom double selection
double_selection(x=X,y=y,d=d, method="hdm")[1]
double_selection(x=X,y=y,d=d, method="gamlr")[1]

