example = function() {
  
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
  lasso_sim(beta.e = beta.e,beta.cx = 5,Kc=50,Ke=5,Ky=50,  Ku=400, n=500)
  
  
    
  # Normal lasso with gamlr works very well but hdm not
  
  lasso_sim(beta.d = 1,beta.e = 1.5,beta.y = 1,Kc=50,Ke=50,Ky=50,  Ku=400, n=500)
  
  
  lasso_sim(beta.d = 1,beta.e = 0.8,beta.y = 1,Kc=5,Ke=5,Ky=5,  Ku=400, n=500)
  
  
  # Where does hdm work better than gamlr?
  lasso_sim(beta.d = 1,beta.e = 0,beta.y = 1.5,sd.d = 10,Kc=30,Ke=1,Ky=500,  Ku=500, n=600)
  
  lasso_sim(beta.d = 1,beta.e = 0,beta.y = 1.5,sd.d = 10,Kc=30,Ke=1,Ky=500,  Ku=500, n=2000)

  
  library(sktools)
  sim = simulation.study(lasso_sim, par=list(beta.d = 1,beta.e = 1.5,beta.y = 1,Kc=50,Ke=50,Ky=50,  Ku=400, n=500), repl=100)  

  # Much more sparsity
  sim3 = simulation.study(lasso_sim, par=list(beta.d = 1,beta.e = 1.5,beta.y = 1,Kc=50,Ke=50,Ky=50,  Ku=1000, n=500), repl=50)  
  
    
  sim2 = simulation.study(lasso_sim, par=list(beta.d = 1,beta.e = 0,beta.y = 1.5,sd.d = 10,Kc=30,Ke=1,Ky=500,  Ku=500, n=600), repl=50)  

  #saveRDS(sim3, "lasso_sim_3.Rds")
  
  #saveRDS(sim2, "lasso_sim_2.Rds")
  
  #saveRDS(sim, "lasso_sim_1.Rds")
  sim = readRDS("lasso_sim_1.Rds")
  sim = readRDS("lasso_sim_2.Rds")
  sim = readRDS("lasso_sim_3.Rds")
  
  sim = sim[,-1]
  library(dplyr)
  sim %>%
    summarize(
      rmse_ols_short = sqrt(mean((ols.short-beta.d)^2)),
      rmse_hdm = sqrt(mean((hdm-beta.d)^2)),
      rmse_gamlr = sqrt(mean((gamlr-beta.d)^2)),
      bias_ols_short = mean((ols.short-beta.d)),
      bias_hdm = mean((hdm-beta.d)),
      bias_gamlr = mean(gamlr-beta.d)
    )
  
  library(ggplot2)
  ggplot(sim, aes(x=beta.gamlr)) + geom_density(fill="red", alpha=0.5)
  
}


lasso_sim = function(
  n = 500,
  K = 500,
  Kc = round(K/4),
  Ke = round(K/4),
  Ku = round(K/4),
  Ky = round(K/4),
  sd.xc = 1,sd.xe  = 1,sd.xu = 1,sd.xy = 1,
  sd.d = 1, sd.y = 1, sd.beta = 0,
  beta.e = 1, beta.cx = 1, beta.cy = 1, beta.y = 1,
  beta.d = 1,
  return.what = c("data","coef","details")[3],
  models = list(
    list(lasso="gamlr",type="free_pl"),
    list(lasso="gamlr",type="pl"),
    list(lasso="rlasso",type="pl"),
    list(lasso="gamlr",type="ds"),
    list(lasso="rlasso",type="ds"),
    list(name="ds_gamma_02", lasso="rlasso",type="ds", args = list(penalty=list(c=1.05,gamma=0.2)))
  )
  
) {
  library(restorepoint); library(dplyr)
  restore.point("lasso_sim")
  
  K.vec = c(Kc=Kc, Ke=Ke,Ku=Ku,Ky=Ky)
  Xc = rnorm(Kc*n, 0, sd.xc) %>% matrix(nrow=n) 
  Xe = rnorm(Ke*n, 0, sd.xe) %>% matrix(nrow=n) 
  Xu = rnorm(Ku*n, 0, sd.xu) %>% matrix(nrow=n) 
  Xy = rnorm(Ky*n, 0, sd.xy) %>% matrix(nrow=n) 
  
  d = Xe %*% rnorm(Ke,beta.e, sd.beta) + 
    Xc %*% rnorm(Kc,beta.cx, sd.beta) + rnorm(n, 0, sd.d)
  y = beta.d*d + Xy %*% rnorm(Ky,beta.y,sd.beta) + Xc %*% rnorm(Kc,beta.cy, sd.beta) + rnorm(n, 0, sd.y)
  eps = y-beta.d*d

  X = cbind(Xc,Xe,Xu,Xy)
  colnames(X) = c(paste0("xc",1:Kc),paste0("xe",1:Ke), paste0("xu",1:Ku), paste0("xy",1:Ky))
  dX = cbind(d,X)
  colnames(dX)[1] = "d"
  
  if (return.what=="data") return(list(y=y,d=d,X=X,dX=dX))
  
  
  library(gamlr)
  library(hdm)
  beta.ols.short = coef(lm.fit(x=cbind(1,d),y=y))[2]
  
  model_names = lapply(models, function(model) {
    if (is.null(model$name)) {
      paste0(model$type,"_", model$lasso)
    } else {
      model$name
    }
  })
  names(models) = model_names
  
  res.li = lapply(models, function(model) {
    estimate_lasso_sim_model(model,y=y,d=d,X=X,dX=dX, details=details)
  })
  
  coefs = unlist(lapply(res.li, function(res) res$coef), use.names=FALSE)
  names(coefs) = names(res.li)
  if (return.what =="coef") {
    add = c(beta.d, beta.ols.short,cor(eps,d))
    names(add) = c("beta.d","ols.short", "cor_eps_d")
    return(c(coefs,add))  
  }
  
  res = res.li[[1]]
  
  counts = bind_rows(lapply(res.li, function(res) {
    vs = vars_counts(res$vars)
    vs1 = vars_counts(res$vars1)
    vs2 = vars_counts(res$vars2)
    names(vs1) = paste0(names(vs1),"1")
    names(vs2) = paste0(names(vs2),"2")
    c(vs,vs1,vs2)
  }))
  
  cbind(quick_df(
    model = names(models),
    coef = coefs,  
    ols.short = beta.ols.short
    ),counts)
}

estimate_lasso_sim_model = function(model, d,y,X,dX, details = TRUE) {
  restore.point("estimate_lasso_sim_model")
  
  # post lasso with d as free variable
  if (model$type == "free_pl") {
    if (model$lasso != "gamlr")
      error("free post lasso estimation is so far only implemented with gamlr")
    args = c(list(x=dX,y=y,free=1), model$args)
    reg = do.call(gamlr, args)
    coefs = post_lasso_coef(reg,dX,y)
    list(coef = coefs[1], vars = names(coefs), vars1=NULL, vars2=names(coefs))
  # post lasso without free variable in intia lasso
  } else if (model$type == "pl") {
    args = c(list(x=dX,y=y), model$args)
    if (model$lasso=="gamlr") {
      reg = do.call(gamlr, args)
    } else if (model$lasso == "rlasso") {
      reg = do.call(rlasso, args)
    }
    coefs = post_lasso_coef(reg,dX,y,add.var = "d")
    list(coef = coefs[1], vars = names(coefs)[-1], vars1=NULL, vars2=names(coefs)[-1])
  } else if (model$type == "ds") {
    args = c(list(d=d,x=X,y=y,lasso.method = model$lasso, just.d.coef = TRUE), model$args)
    res = do.call(double_selection, args)
  } else {
    stop(paste0("Model type ", model$type ," not known. So far only model types ds and free_pl are implemented."))
  }
}

vars_counts = function(vars) {
  if (is.null(vars)) return(c(num.vars=0, xc=0,xe=0,xu=0,xy=0))
  res = c(
    xc = sum(startsWith(vars,"xc")),
    xe=sum(startsWith(vars,"xe")),
    xu=sum(startsWith(vars,"xu")),
    xy = sum(startsWith(vars,"xy"))
  )
  count = sum(res)
  c(num.vars=count,res)
}



vars_shares = function(vars) {
  if (is.null(vars)) return(c(num.vars=0, xc=NA,xe=NA,xu=NA,xy=NA))
  res = c(
    xc = sum(startsWith(vars,"xc")),
    xe=sum(startsWith(vars,"xe")),
    xu=sum(startsWith(vars,"xu")),
    xy = sum(startsWith(vars,"xy"))
  )
  count = sum(res)
  c(num.vars=count,round(res / count,2))

}

