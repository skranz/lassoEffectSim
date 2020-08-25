double.selection = function(d,x,y,..., lasso.method=c("rlasso","gamlr")[1],dvar="d", keep.intercept=FALSE, just.d.coef = FALSE) {
  args = list(...)
  if (lasso.method == "rlasso") {
    library(hdm)
    lasso1 = rlasso(x=x,y=d,...)
    lasso2 = rlasso(x=x,y=y,...)
  } else if (lasso.method=="gamlr") {
    library(gamlr)
    lasso1 = gamlr(x=x,y=d,...)
    lasso2 = gamlr(x=x,y=y,...)
  }
  restore.point("double.selection")
  vars1 = names(lasso_coef(lasso1)) 
  vars2 = names(lasso_coef(lasso2))
  vars = union(vars1,vars2)
  X = cbind(1,d,x[,vars,drop=FALSE])
  if (NCOL(X)>2) {
    colnames(X)[1:2] = c("(Intercept)",dvar)
  } else {
    colnames(X) = c("(Intercept)",dvar)
  }
  co = coef(lm.fit(x=X,y=y))
  if (!keep.intercept) co = co[-1]
  if (just.d.coef) co = co[1]
  list(coef = co,vars=vars, vars1=vars1, vars2=vars2)
}


# Call cv.glmnet with free variables that will not be penalized
cv.glmnet.free = function(x,y,free=NULL, ...) {
  library(glmnet)
  if (is.null(free)) {
    cv.glmnet(x=x,y=y,...)
  } else {
    nvars = NCOL(x)
    penalty.factor = rep(1, nvars)
    if (is.character(free))
      free = match(free, colnames(x))
    penalty.factor[free] = 0
    cv.glmnet(x=x,y=y,penalty.factor = penalty.factor)
  }
}

# Call glmnet with free variables that will not be penalized
glmnet.free = function(x,y,free=NULL, ...) {
  library(glmnet)
  if (is.null(free)) {
    glmnet(x=x,y=y,...)
  } else {
    nvars = NCOL(x)
    penalty.factor = rep(1, nvars)
    if (is.character(free))
      free = match(free, colnames(x))
    penalty.factor[free] = 0
    glmnet(x=x,y=y,penalty.factor = penalty.factor)
  }
}

# Return non-zero lasso coefficients
lasso_coef = function(lasso,..., keep.intercept=FALSE) {
  restore.point("lasso_coef")
  if (is(lasso,"rlasso")) {
    co = coef(lasso)
    if (!keep.intercept) {
      if (names(co)[1]=="(Intercept)") co = co[-1]
      co = co[co!=0]
      return(co)
    }
  } else if (is(lasso, "cv.glmnet")) {
    co = as.matrix(coef(lasso, s = "lambda.min"))
    
  } else if (is(lasso,"glmnet")) {
    stop("Not yet implemented for glmnet. Please call cv.glmnet.")
  } else if (is(lasso, "gamlr")) {
    co = as.matrix(coef(lasso,...))
  }
  rows = co[,1] != 0
  rows[1] = keep.intercept
  rows
  co = co[rows,]
  co  
}

post_lasso_coef = function(lasso, X,y, add.var=NULL, keep.intercept=FALSE) {
  restore.point("post_lasso_coef")
  co = lasso_coef(lasso)
  vars = unique(c(add.var,names(co)))
  X = as.matrix(X[,vars])
  reg = lm.fit(x = cbind(1,X),y=y)
  co = coef(reg)
  if (!keep.intercept) co = co[-1]
  co
}

with.random.seed = function(expr, seed = 1234567890) 
{
  old.seed = get(".Random.seed", .GlobalEnv)
  set.seed(seed)
  ret = eval(expr)
  assign(".Random.seed", old.seed, .GlobalEnv)
  runif(1)
  return(ret)
}

quick.df = function(...) {
  as_tibble(list(...))
}
