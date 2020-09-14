# These functions perform the large simulations
# whose results are shown in the blog

# Both take quite a while to run

control.exo.sim = function() {
  source("lasso_tools.R")
  source("lasso_sim.R")
  set.seed(1)
  
  sim.and.est1 = function() {
    #mat = lasso_sim(alpha=1,n=800,Kc=50,Ke=50,Ky=50,Ku=500,return.what = "data")
    mat = lasso_sim(alpha=1,n=700,Kc=50,Ke=50,Ky=50,Ku=700,return.what = "data")
    
    
    y = mat[,1]
    d = mat[,2]
    X = cbind(1,mat[,-1]); colnames(X)[1] = "const"
    
    xc.cols = 3:52; xe.cols = 53:102
    c(
      coef(lm.fit(y=y,x=X[,c(1:2, xc.cols[1:49])]))[2],
      coef(lm.fit(y=y,x=X[,c(1:2, xc.cols[1:49], xe.cols)]))[2]
    )
  }
  library(dplyr)
  res = replicate(1000,sim.and.est1(),simplify = FALSE) %>% bind_rows()
  colnames(res) = c("alpha.hat.no.xe","alpha.hat.control.xe")
  
  sim = tibble(
    reg=c(rep("dont_control_xe",NROW(res)),rep("control_xe",NROW(res))),
    alpha.hat = c(res[[1]],res[[2]]),
    alpha = 1
  )
  saveRDS(sim,"control_exo_sim.Rds")
  
  sim = readRDS("control_exo_sim.Rds")
  head(sim,3)
  library(ggplot2)
  ggplot(sim, aes(x=alpha.hat, fill=reg)) + geom_density() + facet_wrap(~reg) + geom_vline(xintercept=1, alpha=0.7)  
  
  sim %>%
    group_by(reg) %>%
    summarize(
      bias = mean(alpha.hat-alpha),
      rmse = sqrt(mean((alpha.hat-alpha)^2))  
    )
}


# Takes very long to run
lasso.sel.sim = function() {
  source("lasso_tools.R")
  source("lasso_sim.R")
  
  models1 = lasso_sim_default_models(names=c("gamlr_simple","rlasso_double_sel"))
  
  models2 = list(
    rlasso_double_sel_c106 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1.06))),
    rlasso_double_sel_c100 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1)))    
  )
  
  models = c(models1, models2)
  library(dplyr)
  li = replicate(n=1000, lasso_sim(alpha=1,n=700,Kc=50,Ke=50,Ky=50,Ku=700,return.what = "details",models = models),simplify = FALSE)
  
  sim = bind_rows(li)
  saveRDS(sim, "lasso_sel_sim.Rds")
  
  sim %>%
    group_by(model) %>%
    summarize(
      bias = mean(coef-1),
      se = sd(coef)
    )
  library(ggplot2)
  ggplot(sim, aes(x=coef,fill=model)) + geom_density() +
    facet_wrap(~model, scales="free_y")  
}


# Takes very long to run
lasso.sel.sim = function() {
  source("lasso_tools.R")
  source("lasso_sim.R")
  
  models1 = lasso_sim_default_models(names=c("gamlr_simple","rlasso_double_sel"))
  
  models2 = list(
    rlasso_double_sel_c106 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1.06))),
    rlasso_double_sel_c100 = list(lasso.fun="rlasso",type="double_sel", args=list(penalty=list(c=1)))    
  )
  
  models = c(models1, models2)
  library(dplyr)
  li = replicate(n=1000, lasso_sim(alpha=1,n=100,Kc=10,Ke=10,Ky=5,Ku=20, beta.ed = 10, beta.ey = 0.5, beta.cd = 0.5, beta.cy = 10, models=c("gamlr_simple","rlasso_double_sel"))
                 ,simplify = FALSE)
  
  sim = bind_rows(li)
  saveRDS(sim, "lasso_sim3.Rds")
  
  sim = readRDS("lasso_sim3.Rds")
  sim %>%
    group_by(model) %>%
    summarize(
      bias = mean(coef-1),
      se = sd(coef),
      num.vars = mean(num.vars),
      xe = mean(xe),
      xc = mean(xc)
    )  
  library(ggplot2)
  
  ggplot(sim, aes(x=coef,fill=model)) + geom_density() +
    facet_wrap(~model, scales="free_y")  
}


