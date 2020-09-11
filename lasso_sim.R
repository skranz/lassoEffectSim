example = function() {
  # Try to find zero rlasso selection
  set.seed(9)
  lasso_sim(alpha=1,n=500,Kc=50,Ke=50,Ky=50,Ku=500,sd.noise=0.1)
  
  set.seed(9)
  mat = lasso_sim(alpha=1,n=500,Kc=50,Ke=50,Ky=50,Ku=500,sd.noise=0.1,return.what = "data")
  
  y = mat[,1]
  d = mat[,2]
  X = cbind(1,mat[,-1]); colnames(X)[1] = "const"
  
  lasso1 = rlasso(x=X[,-2],y=d)
  coef(lasso1)[coef(lasso1)!=0]

  lasso2 = rlasso(x=X[,-2],y=y)
  coef(lasso2)[coef(lasso2)!=0]
  
  
  res = double_selection(d=d,y=y,x=X[,-2],lasso.fun = "rlasso")
  res

  # Find a good example without were simple lasso works and
  set.seed(1)
  lasso_sim(alpha=1,n=700,Kc=50,Ke=50,Ky=50,Ku=700,return.what = "details")
  
  n = 700
  penalty = list(homoscedastic = FALSE,
                 X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = 0.1/log(n))
  
  models = list(
    rlasso_double_sel_c106 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1.06))),
    rlasso_double_sel_c100 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1)))    
  )
  set.seed(1)
  lasso_sim(alpha=1,n=700,Kc=50,Ke=50,Ky=50,Ku=700,return.what = "details",models = models)
  rlasso
    
  lasso_sim(alpha=1,n=500,Kc=30,Ke=30,Ky=30,Ku=300,return.what = "details")
  

  cor(mat[,1],mat[,3])
  cor(mat[,1],mat[,4])
  cor(mat[,1],mat[,5])
  cor(mat[,1],mat[,6])
  
  
  beta.e = 2
  # If everything is fairly sparse, it all seems to work
  # including normal gamlr lasso
  lasso_sim(beta.e = beta.e, Kc=5,Ke=5,Ky=5,  Ku=500, n=500)

  # Make uniformely less sparse
  lasso_sim(beta.e = beta.e,Kc=50,Ke=50,Ky=50,  Ku=500, n=500)
  # gamlr still works somewhat but double selection and oml break down
  
  # Just add a lot of exogenous variables
  lasso_sim(beta.e = beta.e,Kc=5,Ke=50,Ky=5,  Ku=400, n=500)
  
  # Just add a lot of confounders: gamlr works worst
  lasso_sim(beta.e = beta.e,Kc=50,Ke=5,Ky=5,  Ku=400, n=500)
  
  # Just add a lot of variables that just affect y
  lasso_sim(beta.e = beta.e,Kc=5,Ke=5,Ky=50,  Ku=400, n=500)
  

  beta.e = 2
  # If everything is fairly sparse, it all seems to work
  # including normal gamlr lasso
  lasso_sim(beta.e = beta.e, Kc=5,Ke=5,Ky=5,  Ku=500, n=500)
  
  # Make uniformely less sparse
  lasso_sim(beta.e = beta.e,Kc=50,Ke=50,Ky=50,  Ku=500, n=500)
  lasso_sim(beta.e = beta.e,Kc=50,Ke=50,Ky=50,  Ku=500, n=500)
  # gamlr still works somewhat but double selection and oml break down
  
  # Just add a lot of exogenous variables
  lasso_sim(beta.e = beta.e,Kc=5,Ke=50,Ky=5,  Ku=400, n=500)

  # Just add a lot of variables that just affect y
  lasso_sim(beta.e = beta.e,Kc=5,Ke=5,Ky=50,  Ku=400, n=500)
  
    
  # Just add a lot of confounders: gamlr works worst if sd.x is small
  lasso_sim(beta.e = beta.e,Kc=50,Ke=5,Ky=5,  Ku=400, n=500)
  lasso_sim(beta.e = beta.e,sd.d=5,Kc=50,Ke=5,Ky=5,  Ku=400, n=500)
  lasso_sim(beta.e = beta.e,sd.d=1,sd.beta=4,Kc=50,Ke=5,Ky=5,  Ku=400, n=500)
  
  # Just add a lot of confounders: gamlr works worst if sd.x is small
  lasso_sim(beta.e = 1,sd.d = 1,Kc=50,Ke=5,Ky=50,  Ku=400, n=500)
  lasso_sim(beta.e = 1,sd.d = 1,sd.beta=4,Kc=50,Ke=5,Ky=50,  Ku=400, n=500)
  lasso_sim(beta.e = 1,sd.d = 10,Kc=50,Ke=1,Ky=50,  Ku=400, n=500)
  
  # Just add a lot of confounders: gamlr works worst
  lasso_sim(beta.e = beta.e,beta.cd = 5,Kc=50,Ke=5,Ky=50,  Ku=400, n=500)
  
  
    
  # Normal lasso with gamlr works very well but hdm not
  
  lasso_sim(alpha = 1,beta.e = 1.5,beta.y = 1,Kc=50,Ke=50,Ky=50,  Ku=400, n=500)
  
  
  lasso_sim(alpha = 1,beta.e = 0.8,beta.y = 1,Kc=5,Ke=5,Ky=5,  Ku=400, n=500)
  
  
  # Where does hdm work better than gamlr?
  lasso_sim(alpha = 1,beta.e = 0,beta.y = 1.5,sd.d = 10,Kc=30,Ke=1,Ky=500,  Ku=500, n=600)
  
  lasso_sim(alpha = 1,beta.e = 0,beta.y = 1.5,sd.d = 10,Kc=30,Ke=1,Ky=500,  Ku=500, n=2000)

  
  # Some simulations where both xe and xc are confounders (last chapter
  # of blog)
  set.seed(1)
  lasso_sim(alpha=1,n=100,Kc=20,Ke=5,Ky=5,Ku=20, beta.e = 20, beta.ey = 0.5, beta.cd = 1, beta.cy = 2)
  
  set.seed(1)
  lasso_sim(alpha=1,n=100,Kc=10,Ke=10,Ky=5,Ku=20, beta.e = 10, beta.ey = 0.5, beta.cd = 0.5, beta.cy = 10)
  
}


lasso_sim = function(
  n = 500,
  K = 500,
  Kc = round(K/4),
  Ke = round(K/4),
  Ku = round(K/4),
  Ky = round(K/4),
  sd.xc = 1,sd.xe  = 1,sd.xu = 1,sd.xy = 1,
  sd.d = 1, sd.y = 1, sd.beta = 0, sd.noise = 0,
  alpha = 1,beta=1,
  beta.e = beta.ed, beta.cd = beta, beta.cy = beta, beta.y = beta, beta.ey = 0, beta.ed=beta,
  return.what = c("data","coef","details")[3],
  models = c("short_ols", "gamlr_simple","rlasso_double_sel","gamlr_double_sel"),
  count.what = "vars"
) {
  library(restorepoint); library(dplyr)
  restore.point("lasso_sim")
  
  K.vec = c(Kc=Kc, Ke=Ke,Ku=Ku,Ky=Ky)
  Xc = rnorm(Kc*n, 0, sd.xc) %>% matrix(nrow=n) 
  Xe = rnorm(Ke*n, 0, sd.xe) %>% matrix(nrow=n) 
  Xu = rnorm(Ku*n, 0, sd.xu) %>% matrix(nrow=n) 
  Xy = rnorm(Ky*n, 0, sd.xy) %>% matrix(nrow=n) 
  
  d = Xe %*% rnorm(Ke,beta.e, sd.beta) + 
    Xc %*% rnorm(Kc,beta.cd, sd.beta) + rnorm(n, 0, sd.d)
  y = alpha*d + Xy %*% rnorm(Ky,beta.y,sd.beta) + Xc %*% rnorm(Kc,beta.cy, sd.beta) + Xe %*% rep(beta.ey,Ke) + rnorm(n, 0, sd.y)
  eps = y-alpha*d

  X = cbind(Xc,Xe,Xu,Xy)
  if (sd.noise > 0) {
    noise = rnorm(length(X),0,sd.noise)
    X = X+noise
  }
  
  colnames(X) = c(paste0("xc",1:Kc),paste0("xe",1:Ke), paste0("xu",1:Ku), paste0("xy",1:Ky))

  if (return.what=="data") {
    mat = cbind(y,d,X)
    colnames(mat)[1:2] = c("y","d")
    return(mat)
  }
  
  lasso_sim_estimate(models=models, y=y,d=d,X=X, count.what=count.what, return.what=return.what)
}



lasso_sim_default_models = function(names) {
  models = list(
    short_ols = list(type="short_ols"),
    gamlr_free_d = list(lasso.fun="gamlr",type="post_lasso", args=list(free=1)),
    gamlr_simple = list(lasso.fun="gamlr",type="post_lasso"),
    rlasso_simple = list(lasso.fun="rlasso",type="pl"),
    rlasso_double_sel = list(lasso.fun="rlasso",type="double_sel"),
    gamlr_double_sel = list(lasso.fun="gamlr",type="double_sel")
  )
  
  if (!missing(names)) {
    unknown = setdiff(names, names(models))
    if (length(unknown) > 0) 
      stop(paste0("No model with name ", unknown[1], " specified in lasso_sim_default_models."))
    models = models[names]
  }
  models
}

lasso_sim_estimate = function(mat, models=lasso_sim_default_models(),y,d,X, return.what ="details", count.what="vars") {
  if (missing(y)) {
    if (missing(mat)) {
      stop("You either have to provide mat or the three arguments y,d,X")
    }
    y = mat[,1]
    d = mat[,2]
    X = mat[, -c(1:2)]
    dX = mat[,-1]
  } else {
    dX = cbind(d,X)
    colnames(dX)[1] = "d"
  }
  restore.point("lasso_sim_estimate")
  
  if (is.character(models)) {
    models = lasso_sim_default_models(models)
  }
  
  
  library(gamlr)
  library(hdm)
  model_names = names(models)
  res.li = lapply(models, function(model) {
    lasso_sim_estimate_model(model,y=y,d=d,X=X,dX=dX)
  })
  
  coefs = unlist(lapply(res.li, function(res) res$coef), use.names=FALSE)
  names(coefs) = names(res.li)
  if (return.what =="coef") {
    return(coefs)  
  }
  
  res = res.li[[1]]
  
  counts = bind_rows(lapply(res.li, function(res) {
    
    vs = if ("vars" %in% count.what) vars_counts(res$vars)
    vs1 = if ("vars1" %in% count.what) vars_counts(res$vars1)
    vs2 = if ("vars2" %in% count.what) vars_counts(res$vars2)
    if (!is.null(vs1)) names(vs1) = paste0(names(vs1),"1")
    if (!is.null(vs2)) names(vs2) = paste0(names(vs2),"2")
    c(vs,vs1,vs2)
  }))
  
  cbind(quick_df(
    model = names(models),
    coef = coefs 
  ),counts)
}

lasso_sim_estimate_model = function(model, d,y,X,dX, details = TRUE) {
  restore.point("estimate_lasso_sim_model")
  
  # post lasso with d as free variable
  # post lasso without free variable in intia lasso
  if (model$type == "post_lasso") {
    args = c(list(x=dX,y=y), model$args)
    if (model$lasso.fun=="gamlr") {
      reg = do.call(gamlr, args)
    } else if (model$lasso.fun == "rlasso") {
      reg = do.call(rlasso, args)
    }
    coefs = post_lasso_coef(reg,dX,y,add.var = "d")
    list(coef = coefs[1], vars = names(coefs)[-1], vars1=NULL, vars2=names(coefs)[-1])
  } else if (model$type == "double_sel") {
    args = c(list(d=d,x=X,y=y,lasso.fun = model$lasso.fun, just.d.coef = TRUE), model$args)
    res = do.call(double_selection, args)
  } else if (model$type == "short_ols") {
    coef = coef(lm.fit(y=y,x=cbind(1,d)))[2]
    list(coef=coef, vars=NULL)             
  } else {
    stop(paste0("Model type ", model$type ," not known. So far only model types double_sel, lasso and post_lasso are implemented."))
  }
}


vars_counts = function(vars) {
  if (is.null(vars)) return(c(num.vars=0, xc=0,xe=0,xu=0,xy=0))
  
  
  if (!is.character(vars)) {
    vars = names(lasso_coef(vars))
  }
  
  res = c(
    xc = sum(startsWith(vars,"xc")),
    xe=sum(startsWith(vars,"xe")),
    xu=sum(startsWith(vars,"xu")),
    xy = sum(startsWith(vars,"xy"))
  )
  count = sum(res)
  c(num.vars=count,res)
}




