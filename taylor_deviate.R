

ps.model.fit = function(ps.ds, fm_fit.s, a){
  # Design matrix
  psa_dat = ps.ds$variables
  R = ps.ds$strata[[names(ps.ds$strata)]]
  ps.w = weights(ps.ds)
  lgtreg = svyglm(fm_fit.s, family = binomial, design = ps.ds)
  design.x = model.matrix(lgtreg)
  gamma_est = summary(lgtreg)$coeff[,1]
  #a = fs
  svy.wt = ps.w[R==0]/a
  # Propensity scores
  p.est = lgtreg$fitted.values
  p.c = p.est[R==1]
  pi.c_est = p.c/(1-p.c)*a
  # Linear propensity scores 
  p_score   = lgtreg$linear.predictors

  # Check Taylor deviate method for gamma estimates
  gamma_w_out = gamma_w(psa_dat = psa_dat, ps.w = ps.w, R = R, p.est = p.est,
                        pi.c_est = pi.c_est, design.x = design.x, svy.wt = svy.wt)
  #gamma_w_out$var_gamma_pps; vcov(lgtreg)
  
  # Save the scores and design matrices for later use
  model.par = list(psa_dat = psa_dat,
                   R = R,
                   p.est = p.est, 
                   p_score = p_score,
                   pi.c_est = pi.c_est,
                   design.x=design.x,
                   svy.wt = svy.wt,
                   ps.w = ps.w,
                   a = a)
  return(list(ps.model = summary(lgtreg), model.par = model.par, gamma_w_out = gamma_w_out))
}

inference_beta = function(ps.model.par, 
                          #fm_fit.lg, 
                          fm_fit.cox, a, 
                          t_star, lambda_star, x0, Lambda_t,
                          post.h, 
                          Nh,
                          d="d",
                          t="t"){
  model.par = ps.model.par$model.par
  psa_dat = model.par$psa_dat
  R = model.par$R
  p.est = model.par$p.est
  p_score = model.par$p_score
  pi.c_est =  model.par$pi.c_est
  design.x = model.par$design.x
  design.x.c = design.x[R==1,]
  design.x.s = design.x[R==0,]
  svy.wt = model.par$svy.wt
  ps.w = model.par$ps.w
  a = model.par$a
  
  samp.c = psa_dat[R==1,]
  p_score.c = p_score[R==1]
  p_score.s = p_score[R==0]
  # IPSW
  samp.c$ipsw = exp(-p_score.c)
  ipsw_gamma = -samp.c$ipsw*design.x.c
  # Estimate lg and cox regression coefficients
  ds.c = svydesign(ids=~1, data=samp.c, weights=~ipsw, nest=TRUE)
  #lg_out = svyglm(fm_fit.lg, family=binomial, design = ds.c)
  #lg_est.ipsw = lg_out$coefficients
  cox_out = svycoxph(fm_fit.cox, design = ds.c)
  x.mtrx = as.matrix(model.matrix(cox_out))
  
  rel_hzd = exp(x.mtrx%*%c(cox_out$coefficients))
  cox_est.ipsw = cox_out$coefficients
  
  # Estimate variance
  #x.mtrx = as.matrix(model.matrix(lg_out))
  p = ncol(x.mtrx)
  
  beta_var_ipsw = beta_var_est(x.mtrx = x.mtrx, 
                               #y = "d", 
                               #mu = lg_out$fitted.values,
                               rel_hzd = rel_hzd,
                               dat = samp.c, 
                               pw_gamma = ipsw_gamma,
                               model.par = model.par,
                               pw = "ipsw"#,
                               #ps.w = ps.w
                               )
  
  #beta_wt.ipsw = beta_var_ipsw$Delta_eta_cox[,1:(p-1)]/ps.w
  ps.w = model.par$ps.w
  beta_wt.ipsw = beta_var_ipsw$Delta_eta_cox[,1:p]/ps.w
  

  lambda_out.ipsw = lambda.ar_w (x.mtrx  = as.matrix(model.matrix(cox_out)), 
                                 dat     = samp.c, 
                                 pw      = "ipsw", 
                                 rel_hzd = rel_hzd, 
                                 beta_wt = beta_wt.ipsw, 
                                 pw_list = list(pw_gamma = ipsw_gamma, 
                                                gamma_wt = ps.model.par$gamma_w_out$gamma_wt)) 
  Delta_lambda.ipsw = ps.w*lambda_out.ipsw$lambda_wt
  Lambda_out.ipsw = Lambda_w(lambda = lambda_out.ipsw$lambda, 
                             lambda_wt = lambda_out.ipsw$lambda_wt, 
                             u = lambda_out.ipsw$u, t_star = t_star)
  Lambda.ipsw = Lambda_out.ipsw$Lambda
  DeltaLambda.ipsw = ps.w*Lambda_out.ipsw$Lambda_wt
  
  var_Lambda.ipsw.pps = diag(cov(DeltaLambda.ipsw[R==1,])*sum(R)+cov(DeltaLambda.ipsw[R==0,])*sum(1-R))
  var_Lambda.ipsw.poisson = diag(t((1-pi.c_est)*DeltaLambda.ipsw[R==1,])%*%DeltaLambda.ipsw[R==1,]+
                                 t((1-1/svy.wt)*DeltaLambda.ipsw[R==0,])%*%DeltaLambda.ipsw[R==0,])
  
  absR_out.ipsw = absR_w(beta_est = cox_est.ipsw, Lambda = Lambda_out.ipsw$Lambda[Lambda_t], x0 =x0, 
                    beta_w = beta_wt.ipsw, Lambda_w = Lambda_out.ipsw$Lambda_wt[,Lambda_t])
  absR.ipsw = absR_out.ipsw$absR
  Delta_absR.ipsw = ps.w*absR_out.ipsw$absR_w
  var_absR.ipsw.poisson = diag(cov(Delta_absR.ipsw[R==1,])*sum(R)+cov(Delta_absR.ipsw[R==0,])*sum(1-R))
  var_absR.ipsw.pps = diag(t((1-pi.c_est)*Delta_absR.ipsw[R==1,])%*%Delta_absR.ipsw[R==1,]+
                           t((1-1/svy.wt)*Delta_absR.ipsw[R==0,])%*%Delta_absR.ipsw[R==0,]) 
  
  LambdaG_out.ipsw = LambdaG_w(ar = lambda_out.ipsw$ar, ar_wt = lambda_out.ipsw$ar_wt, 
                               u = lambda_out.ipsw$u, lambda_star = lambda_star, 
                               t_star = t_star)
  LambdaG.ipsw = LambdaG_out.ipsw$LambdaG
  
  DeltaLambdaG.ipsw =ps.w*LambdaG_out.ipsw$LambdaG_wt
  var_LambdaG.ipsw.pps = diag(cov(DeltaLambdaG.ipsw[R==1,])*sum(R)+cov(DeltaLambdaG.ipsw[R==0,])*sum(1-R))
  var_LambdaG.ipsw.poisson = diag(t((1-pi.c_est)*DeltaLambdaG.ipsw[R==1,])%*%DeltaLambdaG.ipsw[R==1,]+
                                  t((1-1/svy.wt)*DeltaLambdaG.ipsw[R==0,])%*%DeltaLambdaG.ipsw[R==0,])
  
  absRG_out.ipsw = absR_w(beta_est = cox_est.ipsw, Lambda = LambdaG_out.ipsw$LambdaG[Lambda_t], x0 =x0, 
                          beta_w = beta_wt.ipsw, Lambda_w = LambdaG_out.ipsw$LambdaG_wt[,Lambda_t])
  absRG.ipsw = absRG_out.ipsw$absR
  Delta_absRG.ipsw = ps.w*absRG_out.ipsw$absR_w
  var_absRG.ipsw.poisson = diag(cov(Delta_absRG.ipsw[R==1,])*sum(R)+cov(Delta_absRG.ipsw[R==0,])*sum(1-R))
  var_absRG.ipsw.pps = diag(t((1-pi.c_est)*Delta_absRG.ipsw[R==1,])%*%Delta_absRG.ipsw[R==1,]+
                            t((1-1/svy.wt)*Delta_absRG.ipsw[R==0,])%*%Delta_absRG.ipsw[R==0,])
  
  #IPSW+post-stratification
  post = post.f(samp = samp.c[samp.c[,d]==1,], wt0 = "ipsw", post.h =post.h, Nh =Nh*a)
  n_c = nrow(samp.c)
  id_d=c(c(1:n_c)[samp.c[,d]==1], c(1:n_c)[samp.c[,d]==0])
  samp.c$post.ipsw=c(samp.c$ipsw[samp.c[,d]==1]*post$f, samp.c$ipsw[samp.c[,d]==0])[order(id_d)]
  
  ds.c = svydesign(ids=~1, data=samp.c, weights=~post.ipsw)
  cox_out = svycoxph(fm_fit.cox, design = ds.c)
  rel_hzd = exp(x.mtrx%*%c(cox_out$coefficients))
  cox_est.pst.ipsw = cox_out$coefficients
  
  post$f = samp.c$post.ipsw/samp.c$ipsw
  
  beta_var_pst.ipsw = beta_var_est(x.mtrx = x.mtrx, 
                                  # y = "d", 
                                   #mu = lg_out$fitted.values,
                                   rel_hzd = rel_hzd,
                                   dat = samp.c, 
                                   pw_gamma = ipsw_gamma,
                                   model.par = model.par,
                                   pw = "post.ipsw",
                                   #ps.w = ps.w,
                                   post = post)
  beta_wt.pst.ipsw = beta_var_pst.ipsw$Delta_eta_cox[,1:p]/ps.w
  
  lambda_out.pst.ipsw = lambda.ar_w (x.mtrx  = as.matrix(model.matrix(cox_out)), 
                                     dat     = samp.c, 
                                     pw      = "post.ipsw", 
                                     rel_hzd = rel_hzd, 
                                     beta_wt = beta_wt.pst.ipsw, 
                                     pw_list = list(pw_gamma = ipsw_gamma, 
                                                    gamma_wt = ps.model.par$gamma_w_out$gamma_wt),
                                     post=post) 
  Lambda_out.pst.ipsw = Lambda_w(lambda = lambda_out.pst.ipsw$lambda, 
                                 lambda_wt = lambda_out.pst.ipsw$lambda_wt, 
                                 u = lambda_out.pst.ipsw$u, t_star = t_star)
  Lambda.pst.ipsw = Lambda_out.pst.ipsw$Lambda
  DeltaLambda.pst.ipsw = ps.w*Lambda_out.pst.ipsw$Lambda_wt
  
  var_Lambda.pst.ipsw.pps = diag(cov(DeltaLambda.pst.ipsw[R==1,])*sum(R)+cov(DeltaLambda.pst.ipsw[R==0,])*sum(1-R))
  var_Lambda.pst.ipsw.poisson = diag(t((1-pi.c_est)*DeltaLambda.pst.ipsw[R==1,])%*%DeltaLambda.pst.ipsw[R==1,]+
                                     t((1-1/svy.wt)*DeltaLambda.pst.ipsw[R==0,])%*%DeltaLambda.pst.ipsw[R==0,])
  
  absR_out.pst.ipsw = absR_w(beta_est = cox_est.pst.ipsw, Lambda = Lambda_out.pst.ipsw$Lambda[Lambda_t], x0 =x0, 
                         beta_w = beta_wt.pst.ipsw, Lambda_w = Lambda_out.pst.ipsw$Lambda_wt[,Lambda_t])
  absR.pst.ipsw = absR_out.pst.ipsw$absR
  Delta_absR.pst.ipsw = ps.w*absR_out.pst.ipsw$absR_w
  var_absR.pst.ipsw.poisson = diag(cov(Delta_absR.pst.ipsw[R==1,])*sum(R)+cov(Delta_absR.pst.ipsw[R==0,])*sum(1-R))
  var_absR.pst.ipsw.pps = diag(t((1-pi.c_est)*Delta_absR.pst.ipsw[R==1,])%*%Delta_absR.pst.ipsw[R==1,]+
                               t((1-1/svy.wt)*Delta_absR.pst.ipsw[R==0,])%*%Delta_absR.pst.ipsw[R==0,]) 
  
  LambdaG_out.pst.ipsw = LambdaG_w(ar = lambda_out.pst.ipsw$ar, ar_wt = lambda_out.pst.ipsw$ar_wt, 
                               u = lambda_out.pst.ipsw$u, lambda_star = lambda_star, 
                               t_star = t_star)
  LambdaG.pst.ipsw = LambdaG_out.pst.ipsw$LambdaG
  
  DeltaLambdaG.pst.ipsw =ps.w*LambdaG_out.pst.ipsw$LambdaG_wt
  var_LambdaG.pst.ipsw.pps = diag(cov(DeltaLambdaG.pst.ipsw[R==1,])*sum(R)+cov(DeltaLambdaG.pst.ipsw[R==0,])*sum(1-R))
  var_LambdaG.pst.ipsw.poisson = diag(t((1-pi.c_est)*DeltaLambdaG.pst.ipsw[R==1,])%*%DeltaLambdaG.pst.ipsw[R==1,]+
                                      t((1-1/svy.wt)*DeltaLambdaG.pst.ipsw[R==0,])%*%DeltaLambdaG.pst.ipsw[R==0,])
  
  absRG_out.pst.ipsw = absR_w(beta_est = cox_est.pst.ipsw, Lambda = LambdaG_out.pst.ipsw$LambdaG[Lambda_t], x0 =x0, 
                          beta_w = beta_wt.pst.ipsw, Lambda_w = LambdaG_out.pst.ipsw$LambdaG_wt[,Lambda_t])
  absRG.pst.ipsw = absRG_out.pst.ipsw$absR
  Delta_absRG.pst.ipsw = ps.w*absRG_out.pst.ipsw$absR_w
  var_absRG.pst.ipsw.poisson = diag(cov(Delta_absRG.pst.ipsw[R==1,])*sum(R)+cov(Delta_absRG.pst.ipsw[R==0,])*sum(1-R))
  var_absRG.pst.ipsw.pps = diag(t((1-pi.c_est)*Delta_absRG.pst.ipsw[R==1,])%*%Delta_absRG.pst.ipsw[R==1,]+
                              t((1-1/svy.wt)*Delta_absRG.pst.ipsw[R==0,])%*%Delta_absRG.pst.ipsw[R==0,])
  

  # KW
  sgn_dist_mtx = outer(p_score.s, p_score.c, FUN = "-")
  # Bandwidths and weights and  # Kernal smoothing adjustment weights
  h = bw.nrd0(p_score.c)
  out = krnwt(sgn_dist_mtx, h, svy.wt = svy.wt, krnfun=dnorm, rm.s = F,
              w_beta=T, design.x.c=design.x.c, design.x.s=design.x.s)
  samp.c$kw = out$psd.wt
  kw_gamma = out$pw_beta
  colnames(kw_gamma) = colnames(ipsw_gamma)
  # Estimate lg and cox regression coefficients
  ds.c = svydesign(ids=~1, data=samp.c, weights=~kw, nest=TRUE)
  #lg_out = svyglm(fm_fit.lg, family=binomial, design = ds.c)
  #lg_est.kw = lg_out$coefficients
  cox_out = svycoxph(fm_fit.cox, design = ds.c)
  rel_hzd = exp(x.mtrx%*%c(cox_out$coefficients))
  cox_est.kw = cox_out$coefficients
  # Estimate variance
  beta_var_kw = beta_var_est(x.mtrx = x.mtrx, 
                            # y = "d", 
                             #mu = lg_out$fitted.values,
                             rel_hzd = rel_hzd,
                             dat = samp.c, 
                             pw_gamma = kw_gamma,
                             model.par = model.par,
                             pw = "kw"#,
                             #ps.w = ps.w
                             )
  #beta_wt.kw = beta_var_kw$Delta_eta_cox[,1:(p-1)]/ps.w
  beta_wt.kw = beta_var_kw$Delta_eta_cox[,1:p]/ps.w
  
  lambda_out.kw = lambda.ar_w (x.mtrx  = as.matrix(model.matrix(cox_out)), 
                                 dat     = samp.c, 
                                 pw      = "kw", 
                                 rel_hzd = rel_hzd, 
                                 beta_wt = beta_wt.kw, 
                                 pw_list = list(pw_gamma = kw_gamma, 
                                                gamma_wt = ps.model.par$gamma_w_out$gamma_wt)
  ) 
  Delta_lambda.kw = ps.w*lambda_out.kw$lambda_wt
  Lambda_out.kw = Lambda_w(lambda = lambda_out.kw$lambda, 
                             lambda_wt = lambda_out.kw$lambda_wt, 
                             u = lambda_out.kw$u, t_star = t_star)
  Lambda.kw = Lambda_out.kw$Lambda
  DeltaLambda.kw = ps.w*Lambda_out.kw$Lambda_wt
  
  var_Lambda.kw.pps = diag(cov(DeltaLambda.kw[R==1,])*sum(R)+cov(DeltaLambda.kw[R==0,])*sum(1-R))
  var_Lambda.kw.poisson = diag(t((1-pi.c_est)*DeltaLambda.kw[R==1,])%*%DeltaLambda.kw[R==1,]+
                               t((1-1/svy.wt)*DeltaLambda.kw[R==0,])%*%DeltaLambda.kw[R==0,])
  
  absR_out.kw = absR_w(beta_est = cox_est.kw, Lambda = Lambda_out.kw$Lambda[Lambda_t], x0 =x0, 
                         beta_w = beta_wt.kw, Lambda_w = Lambda_out.kw$Lambda_wt[,Lambda_t])
  absR.kw = absR_out.kw$absR
  Delta_absR.kw = ps.w*absR_out.kw$absR_w
  var_absR.kw.poisson = diag(cov(Delta_absR.kw[R==1,])*sum(R)+cov(Delta_absR.kw[R==0,])*sum(1-R))
  var_absR.kw.pps = diag(t((1-pi.c_est)*Delta_absR.kw[R==1,])%*%Delta_absR.kw[R==1,]+
                         t((1-1/svy.wt)*Delta_absR.kw[R==0,])%*%Delta_absR.kw[R==0,]) 
  
  LambdaG_out.kw = LambdaG_w(ar = lambda_out.kw$ar, ar_wt = lambda_out.kw$ar_wt, 
                             u = lambda_out.kw$u, lambda_star = lambda_star, 
                             t_star = t_star)
  LambdaG.kw = LambdaG_out.kw$LambdaG
  
  DeltaLambdaG.kw =ps.w*LambdaG_out.kw$LambdaG_wt
  var_LambdaG.kw.pps = diag(cov(DeltaLambdaG.kw[R==1,])*sum(R)+cov(DeltaLambdaG.kw[R==0,])*sum(1-R))
  var_LambdaG.kw.poisson = diag(t((1-pi.c_est)*DeltaLambdaG.kw[R==1,])%*%DeltaLambdaG.kw[R==1,]+
                                t((1-1/svy.wt)*DeltaLambdaG.kw[R==0,])%*%DeltaLambdaG.kw[R==0,])
  
  absRG_out.kw = absR_w(beta_est = cox_est.kw, Lambda = LambdaG_out.kw$LambdaG[Lambda_t], x0 =x0, 
                        beta_w = beta_wt.kw, Lambda_w = LambdaG_out.kw$LambdaG_wt[,Lambda_t])
  absRG.kw = absRG_out.kw$absR
  Delta_absRG.kw = ps.w*absRG_out.kw$absR_w
  var_absRG.kw.poisson = diag(cov(Delta_absRG.kw[R==1,])*sum(R)+cov(Delta_absRG.kw[R==0,])*sum(1-R))
  var_absRG.kw.pps = diag(t((1-pi.c_est)*Delta_absRG.kw[R==1,])%*%Delta_absRG.kw[R==1,]+
                          t((1-1/svy.wt)*Delta_absRG.kw[R==0,])%*%Delta_absRG.kw[R==0,]) 

  #KW+post-stratification
  post = post.f(samp = samp.c[samp.c[,d]==1,], wt0 = "kw", post.h =post.h, Nh =Nh)
  #n_c = nrow(samp.c)
  #id_d=c(c(1:n_c)[samp.c$d==1], c(1:n_c)[samp.c$d==0])
  samp.c$post.kw=c(samp.c$kw[samp.c[,d]==1]*post$f, samp.c$kw[samp.c[,d]==0])[order(id_d)]
  
  ds.c = svydesign(ids=~1, data=samp.c, weights=~post.kw)
  cox_out = svycoxph(fm_fit.cox, design = ds.c)
  rel_hzd = exp(model.matrix(cox_out)%*%c(cox_out$coefficients))
  cox_est.pst.kw = cox_out$coefficients
  
  post$f = samp.c$post.kw/samp.c$kw
  
  beta_var_pst.kw = beta_var_est(x.mtrx = x.mtrx, 
                                   #y = "d", 
                                   #mu = lg_out$fitted.values,
                                   rel_hzd = rel_hzd,
                                   dat = samp.c, 
                                   pw_gamma = kw_gamma,
                                   model.par = model.par,
                                   pw = "post.kw",
                                   #ps.w = ps.w,
                                   post = post)
  beta_wt.pst.kw = beta_var_pst.kw$Delta_eta_cox[,1:p]/ps.w
  
  lambda_out.pst.kw = lambda.ar_w (x.mtrx  = as.matrix(model.matrix(cox_out)), 
                                     dat     = samp.c, 
                                     pw      = "post.kw", 
                                     rel_hzd = rel_hzd, 
                                     beta_wt = beta_wt.pst.kw, 
                                     pw_list = list(pw_gamma = kw_gamma, 
                                                    gamma_wt = ps.model.par$gamma_w_out$gamma_wt),
                                     post=post) 
  Delta_lambda.pst.kw = ps.w*lambda_out.pst.kw$lambda_wt
  Lambda_out.pst.kw = Lambda_w(lambda = lambda_out.pst.kw$lambda, 
                                 lambda_wt = lambda_out.pst.kw$lambda_wt, 
                                 u = lambda_out.pst.kw$u, t_star = t_star)
  Lambda.pst.kw = Lambda_out.pst.kw$Lambda
  DeltaLambda.pst.kw = ps.w*Lambda_out.pst.kw$Lambda_wt
  
  var_Lambda.pst.kw.pps = diag(cov(DeltaLambda.pst.kw[R==1,])*sum(R)+cov(DeltaLambda.pst.kw[R==0,])*sum(1-R))
  var_Lambda.pst.kw.poisson = diag(t((1-pi.c_est)*DeltaLambda.pst.kw[R==1,])%*%DeltaLambda.pst.kw[R==1,]+
                                       t((1-1/svy.wt)*DeltaLambda.pst.kw[R==0,])%*%DeltaLambda.pst.kw[R==0,])
  
  absR_out.pst.kw = absR_w(beta_est = cox_est.pst.kw, Lambda = Lambda_out.pst.kw$Lambda[Lambda_t], x0 =x0, 
                             beta_w = beta_wt.pst.kw, Lambda_w = Lambda_out.pst.kw$Lambda_wt[,Lambda_t])
  absR.pst.kw = absR_out.pst.kw$absR
  Delta_absR.pst.kw = ps.w*absR_out.pst.kw$absR_w
  var_absR.pst.kw.poisson = diag(cov(Delta_absR.pst.kw[R==1,])*sum(R)+cov(Delta_absR.pst.kw[R==0,])*sum(1-R))
  var_absR.pst.kw.pps = diag(t((1-pi.c_est)*Delta_absR.pst.kw[R==1,])%*%Delta_absR.pst.kw[R==1,]+
                                 t((1-1/svy.wt)*Delta_absR.pst.kw[R==0,])%*%Delta_absR.pst.kw[R==0,]) 
  
  LambdaG_out.pst.kw = LambdaG_w(ar = lambda_out.pst.kw$ar, ar_wt = lambda_out.pst.kw$ar_wt, 
                                   u = lambda_out.pst.kw$u, lambda_star = lambda_star, 
                                   t_star = t_star)
  LambdaG.pst.kw = LambdaG_out.pst.kw$LambdaG
  
  DeltaLambdaG.pst.kw =ps.w*LambdaG_out.pst.kw$LambdaG_wt
  var_LambdaG.pst.kw.pps = diag(cov(DeltaLambdaG.pst.kw[R==1,])*sum(R)+cov(DeltaLambdaG.pst.kw[R==0,])*sum(1-R))
  var_LambdaG.pst.kw.poisson = diag(t((1-pi.c_est)*DeltaLambdaG.pst.kw[R==1,])%*%DeltaLambdaG.pst.kw[R==1,]+
                                        t((1-1/svy.wt)*DeltaLambdaG.pst.kw[R==0,])%*%DeltaLambdaG.pst.kw[R==0,])
  
  absRG_out.pst.kw = absR_w(beta_est = cox_est.pst.kw, Lambda = LambdaG_out.pst.kw$LambdaG[Lambda_t], x0 =x0, 
                              beta_w = beta_wt.pst.kw, Lambda_w = LambdaG_out.pst.kw$LambdaG_wt[,Lambda_t])
  absRG.pst.kw = absRG_out.pst.kw$absR
  Delta_absRG.pst.kw = ps.w*absRG_out.pst.kw$absR_w
  var_absRG.pst.kw.poisson = diag(cov(Delta_absRG.pst.kw[R==1,])*sum(R)+cov(Delta_absRG.pst.kw[R==0,])*sum(1-R))
  var_absRG.pst.kw.pps = diag(t((1-pi.c_est)*Delta_absRG.pst.kw[R==1,])%*%Delta_absRG.pst.kw[R==1,]+
                                  t((1-1/svy.wt)*Delta_absRG.pst.kw[R==0,])%*%Delta_absRG.pst.kw[R==0,])
  
  #lg_var.ipsw = rbind(diag(beta_var_ipsw$beta_lg_var.pps)[1:p],
  #                    diag(beta_var_ipsw$beta_lg_var.poisson)[1:p],
  #                    diag(beta_var_ipsw$beta_lg_2stpvar.pps),
  #                    diag(beta_var_ipsw$beta_lg_2stpvar.poisson))
  #lg_var.kw   = rbind(diag(beta_var_kw$beta_lg_var.pps)[1:p],
  #                    diag(beta_var_kw$beta_lg_var.poisson)[1:p],
  #                    diag(beta_var_kw$beta_lg_2stpvar.pps),
  #                    diag(beta_var_kw$beta_lg_2stpvar.poisson))
  #cox_var.ipsw = rbind(diag(beta_var_ipsw$beta_cox_var.pps)[1:(p-1)],
  #                     diag(beta_var_ipsw$beta_cox_var.poisson)[1:(p-1)],
  #                     diag(beta_var_ipsw$beta_cox_2stpvar.pps),
  #                     diag(beta_var_ipsw$beta_cox_2stpvar.poisson))
  #cox_var.kw   = rbind(diag(beta_var_kw$beta_cox_var.pps)[1:(p-1)],
  #                     diag(beta_var_kw$beta_cox_var.poisson)[1:(p-1)],
  #                     diag(beta_var_kw$beta_cox_2stpvar.pps),
  #                     diag(beta_var_kw$beta_cox_2stpvar.poisson))
  cox_var.ipsw = rbind(diag(beta_var_ipsw$beta_cox_var.pps)[1:p],
                       diag(beta_var_ipsw$beta_cox_var.poisson)[1:p])
  cox_var.kw   = rbind(diag(beta_var_kw$beta_cox_var.pps)[1:p],
                       diag(beta_var_kw$beta_cox_var.poisson)[1:p])
  #row.names(lg_var.ipsw)  = paste(c("pps", "poisson", "pps2", "poisson2"), "ipsw", sep="_")
  #row.names(cox_var.ipsw) = paste(c("pps", "poisson", "pps2", "poisson2"), "ipsw", sep="_")
  #row.names(lg_var.kw)  = paste(c("pps", "poisson", "pps2", "poisson2"), "kw", sep="_")
  #row.names(cox_var.kw) = paste(c("pps", "poisson", "pps2", "poisson2"), "kw", sep="_")
  row.names(cox_var.ipsw) = paste(c("pps", "poisson"), "ipsw", sep="_")
  row.names(cox_var.kw)   = paste(c("pps", "poisson"), "kw", sep="_")
  cox_var.pst.ipsw = rbind(diag(beta_var_pst.ipsw$beta_cox_var.pps)[1:p],
                           diag(beta_var_pst.ipsw$beta_cox_var.poisson)[1:p])
  cox_var.pst.kw   = rbind(diag(beta_var_pst.kw$beta_cox_var.pps)[1:p],
                           diag(beta_var_pst.kw$beta_cox_var.poisson)[1:p])
  row.names(cox_var.pst.ipsw) = paste(c("pps", "poisson"), "pst.ipsw", sep="_")
  row.names(cox_var.pst.kw)   = paste(c("pps", "poisson"), "pst.kw", sep="_")
  

  var.Lambda.pps      = rbind(var_Lambda.ipsw.pps,          var_Lambda.kw.pps,
                              var_Lambda.pst.ipsw.pps,      var_Lambda.pst.kw.pps)
  var.Lambda.poisson  = rbind(var_Lambda.ipsw.poisson,      var_Lambda.kw.poisson,
                              var_Lambda.pst.ipsw.poisson,  var_Lambda.pst.kw.poisson)
  var.LambdaG.pps     = rbind(var_LambdaG.ipsw.pps,         var_LambdaG.kw.pps,
                              var_LambdaG.pst.ipsw.pps,     var_LambdaG.pst.kw.pps)
  var.LambdaG.poisson = rbind(var_LambdaG.ipsw.poisson,     var_LambdaG.kw.poisson,
                              var_LambdaG.pst.ipsw.poisson, var_LambdaG.pst.kw.poisson)

  
  est.Lambda  = rbind(Lambda.ipsw,  Lambda.kw,  Lambda.pst.ipsw,  Lambda.pst.kw)
  est.LambdaG = rbind(LambdaG.ipsw, LambdaG.kw, LambdaG.pst.ipsw, LambdaG.pst.kw)
  colnames(est.Lambda) = colnames(est.LambdaG) = colnames(var.Lambda.pps)
  row.names(var.Lambda.pps) = row.names(var.Lambda.poisson) = row.names(est.Lambda) = c("ipsw", "kw", "pst.ipsw", "pst.kw")
  row.names(var.LambdaG.pps) = row.names(var.LambdaG.poisson) = row.names(est.LambdaG) = c("ipsw", "kw", "pst.ipsw", "pst.kw")
  
  est.absR  = rbind(t(absR.ipsw),  t(absR.kw),  t(absR.pst.ipsw),  t(absR.pst.kw))
  est.absRG = rbind(t(absRG.ipsw), t(absRG.kw), t(absRG.pst.ipsw), t(absRG.pst.kw))
  var.absR.pps      = rbind(var_absR.ipsw.pps,          var_absR.kw.pps,
                            var_absR.pst.ipsw.pps,      var_absR.pst.kw.pps)
  var.absR.poisson  = rbind(var_absR.ipsw.poisson,      var_absR.kw.poisson,
                            var_absR.pst.ipsw.poisson,  var_absR.pst.kw.poisson)
  var.absRG.pps     = rbind(var_absRG.ipsw.pps,         var_absRG.kw.pps,
                            var_absRG.pst.ipsw.pps,     var_absRG.pst.kw.pps)
  var.absRG.poisson = rbind(var_absRG.ipsw.poisson,     var_absRG.kw.poisson,
                            var_absRG.pst.ipsw.poisson, var_absRG.pst.kw.poisson)
  row.names(var.absR.pps)  = row.names(var.absR.poisson)  = row.names(est.absR)  = c("ipsw", "kw", "pst.ipsw", "pst.kw")
  row.names(var.absRG.pps) = row.names(var.absRG.poisson) = row.names(est.absRG) = c("ipsw", "kw", "pst.ipsw", "pst.kw")
  
  return(list(#lg_est       = rbind(lg_est.ipsw, lg_est.kw),
              #lg_var.ipsw  = lg_var.ipsw,
              #lg_var.kw    = lg_var.kw,
              cox_est      = rbind(cox_est.ipsw, cox_est.kw, cox_est.pst.ipsw, cox_est.pst.kw),
              cox_var      = rbind(cox_var.ipsw, cox_var.kw, cox_var.pst.ipsw, cox_var.pst.kw),
              ipsw_gamma   = ipsw_gamma,
              kw_gamma     = kw_gamma,
              beta_wt.ipsw = beta_wt.ipsw,
              beta_wt.kw   = beta_wt.kw,
              beta_wt.pst.ipsw = beta_wt.pst.ipsw,
              beta_wt.pst.kw   = beta_wt.pst.kw,
              est.Lambda   = est.Lambda,
              est.LambdaG  = est.LambdaG,
              var.Lambda.pps       = var.Lambda.pps,
              var.Lambda.poisson   = var.Lambda.poisson,
              var.LambdaG.pps      = var.LambdaG.pps,
              var.LambdaG.poisson  = var.LambdaG.poisson,
              est.absR   = est.absR,
              est.absRG  = est.absRG,
              var.absR.pps       = var.absR.pps,
              var.absR.poisson   = var.absR.poisson,
              var.absRG.pps      = var.absRG.pps,
              var.absRG.poisson  = var.absRG.poisson))
}


#beta_w=residuals(cox_out, "dfbeta", weighted=TRUE)

beta_var_est = function(d="d",
                        t="t",
                        x.mtrx , 
                        rel_hzd ,
                        dat, 
                        pw_gamma ,
                        model.par,
                        pw,
                        post=NULL
                        ){
  a = model.par$a
  psa_dat  = model.par$psa_dat
  R        = model.par$R
  p.est    = model.par$p.est
  p.c      = p.est[R==1]
  p.s      = p.est[R==0]
  design.x = model.par$design.x
  q = ncol(design.x)
  design.x.c = design.x[R==1,]
  design.x.s = design.x[R==0,]
  svy.wt = model.par$svy.wt
  n_s = sum(R==0)
  pi.c_est = model.par$pi.c_est
  p = ncol(x.mtrx)
  ps.w = model.par$ps.w
  #Derivative of eta=(beta, gamma) w.r.t. the original weight (1 for cohort, w for survey sample)
  eta_w = function(U_list){
    Ui_pw  = U_list$Ui_pw
    U_beta = U_list$U_beta
    x.mtrx = U_list$x.mtrx
    S_gamma = -t(ps.w*p.est*(1-p.est)*design.x)%*%design.x
    U_gamma = t(Ui_pw)%*%pw_gamma
    U_beta_inv = solve(U_beta)
    S_gamma_inv = solve(S_gamma)
    b=-U_beta_inv%*%U_gamma%*%S_gamma_inv
    S_beta = matrix(0, nrow=ncol(design.x.c), ncol=ncol(x.mtrx))
    phi_inv = rbind(cbind(U_beta_inv, b),
                    cbind(S_beta, S_gamma_inv))
    
    Si = (R-p.est)*design.x

    eta_w = -t(phi_inv%*%t(cbind(rbind(dat[,pw]*Ui_pw, matrix(0, n_s, ncol(x.mtrx))),
                                Si)))
    
    Delta_eta = ps.w*eta_w
    var_eta_pps = cov(Delta_eta[R==1,])*sum(R)+cov(Delta_eta[R==0,])*sum(1-R)
    var_eta_poisson = t((1-pi.c_est)*Delta_eta[R==1,])%*%Delta_eta[R==1,]+
      t((1-1/svy.wt)*Delta_eta[R==0,])%*%Delta_eta[R==0,]
    return(list(Delta_eta = Delta_eta,
                var_eta_pps = var_eta_pps,
                var_eta_poisson = var_eta_poisson))
  }
  
  # two-step variance
  #two_step.var = function(var_beta.pw, beta_pw, var_gamma){
  #  v2stp.1 = var_beta.pw  # assume SRS for cohort sampling
  #  beta_gamma = t(beta_pw)%*%pw_gamma
  #  v2stp.2 =  beta_gamma%*%var_gamma%*%t(beta_gamma)
  #  v2.stp = v2stp.1+v2stp.2
  #  v2.stp
  #}
  
  #beta_pw_lg_out = beta_pw.lg(x.mtrx = x.mtrx, 
  #                            mu = mu, 
  #                            dat = dat, 
  #                            y = y, 
  #                            pw = pw, 
  #                            pi.c_est = pi.c_est)
  #beta_pw_lg_out$var_beta_pps.pw; vcov(lg_out)
  #eta_w_lg_out = eta_w(beta_pw_lg_out$U_list)
  #beta_lg_2stp_pps = two_step.var(var_beta.pw = beta_pw_lg_out$var_beta_pps.pw, 
  #                                beta_pw = beta_pw_lg_out$beta_pw, 
  #                                var_gamma = eta_w_lg_out$var_eta_pps[-(1:p), -(1:p)])
  #beta_lg_2stp_poisson = two_step.var(var_beta.pw = beta_pw_lg_out$var_beta_poisson.pw, 
  #                                beta_pw = beta_pw_lg_out$beta_pw, 
  #                                var_gamma = eta_w_lg_out$var_eta_poisson[-(1:p), -(1:p)])
  
  beta_pw_cox_out = beta_pw.cox(d=d,
                                t=t,
                                #x.mtrx = as.matrix(x.mtrx[,-1]), 
                                x.mtrx = as.matrix(x.mtrx), 
                                rel_hzd = rel_hzd, 
                                dat = dat, 
                                pw = pw, 
                                pi.c_est = pi.c_est,
                                post=post)
  #beta_pw_cox_out$var_beta_pps.pw; #cox_out$var
  eta_w_cox_out = eta_w(beta_pw_cox_out$U_list)
  #beta_cox_2stp_pps = two_step.var(var_beta.pw = beta_pw_cox_out$var_beta_pps.pw, 
  #                                 beta_pw = beta_pw_cox_out$beta_pw, 
  #                                 var_gamma = eta_w_cox_out$var_eta_pps[-(1:(p-1)), -(1:(p-1))])
  #beta_cox_2stp_poisson = two_step.var(var_beta.pw = beta_pw_cox_out$var_beta_poisson.pw, 
  #                                     beta_pw = beta_pw_cox_out$beta_pw, 
  #                                     var_gamma = eta_w_cox_out$var_eta_poisson[-(1:(p-1)), -(1:(p-1))])
  
  return(list(Delta_eta_cox            = eta_w_cox_out$Delta_eta,
              beta_cox_var.pps         = eta_w_cox_out$var_eta_pps,
              beta_cox_var.poisson     = eta_w_cox_out$var_eta_poisson#,
              #beta_cox_2stpvar.pps     = beta_cox_2stp_pps ,
              #beta_cox_2stpvar.poisson = beta_cox_2stp_poisson,
              #beta_lg_var.pps          = eta_w_lg_out$var_eta_pps,
              #beta_lg_var.poisson      = eta_w_lg_out$var_eta_poisson,
              #beta_lg_2stpvar.pps      = beta_lg_2stp_pps ,
              #beta_lg_2stpvar.poisson  = beta_lg_2stp_poisson
              ))
}


######################################################################

# Derivative of S and gamma (propensity model) w.r.t. original weight
gamma_w = function(psa_dat, ps.w, R, p.est, pi.c_est, design.x, svy.wt){
  Si = (R-p.est)*design.x
  S_gamma = -t(ps.w*p.est*(1-p.est)*design.x)%*%design.x
  gamma_wt = -Si%*%solve(S_gamma)
  Delta_gamma =  ps.w*gamma_wt
  
  var_gamma_pps = cov(Delta_gamma[R==1,])*sum(R)+cov(Delta_gamma[R==0,])*sum(1-R)
  var_gamma_poisson = t((1-pi.c_est)*Delta_gamma[R==1,])%*%Delta_gamma[R==1,]+
    t((1-1/svy.wt)*Delta_gamma[R==0,])%*%Delta_gamma[R==0,]
  return(list(gamma_wt = gamma_wt,
              Delta_gamma = Delta_gamma,
              var_gamma_poisson = var_gamma_poisson,
              var_gamma_pps = var_gamma_pps)
  )
}

################################ FUNCTION f_w_mtrx ##################################################
# f_w_mtrx is a function 
f_w_mtrx = function(f_w, mtrx){
  mtrx = as.matrix(mtrx)
  if((ncol(f_w)==nrow(mtrx))){
    out = as.matrix(f_w)%*%mtrx
  }else{
    mtrx.o = as.matrix(mtrx[order(f_w$post.var), ])
    id.o  = (1:nrow(mtrx))[order(f_w$post.var)]
    f_w.o =aggregate(f_w[,1], list(f_w[,2]), mean)[,2]
    n.grp = length(f_w.o)
    freq_post.var = table(f_w$post.var)
    up = cumsum(freq_post.var)
    lw = c(0, up[-n.grp])+1
    out=NULL
    for(i in 1:n.grp){
      #out = rbind(out, matrix(f_w.o[i], freq_post.var[i], freq_post.var[i])%*%mtrx.o[lw[i]:up[i],])
      out = rbind(out, 
                  t(matrix(as.matrix(f_w.o[i]*colSums(matrix(mtrx.o[lw[i]:up[i],],freq_post.var[i],)))[,rep(1,freq_post.var[i])],
                           ,freq_post.var[i],byrow=F)))
      #print(i)
    }
    out=as.matrix(out[order(id.o),])
  }
  out
}


#Derivative of U and beta (cox regression) w.r.t. pseudo weight
beta_pw.cox = function(x.mtrx, rel_hzd, dat, pw, t="t", d="d", pi.c_est=NULL, post=NULL){
  if(!is.null(post)){
    dat$f = post$f
    f_w = post$f_w
    dat[,pw] = dat[,pw]/dat$f
  }else{dat$f=1}
  n = nrow(dat)
  p = ncol(x.mtrx)
  if(is.null(pi.c_est)) pi.c_est = 1/dat[,pw]
  dat$id = 1:n
  dat$rel_hzd = rel_hzd
  dat$pw_e = dat$f*dat[,pw]*dat$rel_hzd
  
  x.mtrx = as.matrix(x.mtrx[order(dat[,t]),])
  dat = dat[order(dat[,t]),]
  dat$H_dnom = rev(cumsum(rev(dat$pw_e)))# H2
  dat$H_num = as.matrix(cumsum(as.data.frame(c(dat$pw_e)*x.mtrx)[n:1,]))[n:1,] #H1
  
  ties = duplicated(dat[,t])
  if(sum(ties)>0){
    H_uniq = dat[!ties,c(t, "H_dnom", "H_num")]
    H_dat = H_uniq[rep(1:nrow(H_uniq), table(dat[,t])),]
    #h_dat[,t]==dat[,t]
    dat[, c("H_dnom", "H_num")] = H_dat[, c("H_dnom", "H_num")]
    
  }
  
  x.mtrx = as.matrix(x.mtrx[order(dat$id), ])
  dat = dat[order(dat$id), ]
  
  H = as.matrix(dat$H_num/dat$H_dnom)
  
  d_indx = which(dat[,d]==1)
  dat1 = dat[d_indx,]
  # create a two-column matrix sum_wt.t, with the first column being the unique event time, 
  # the second column being the sum of weights for individuals having the event time t
  if(sum(duplicated(dat1[, t]))>0){
    dat1$weight = dat1[,pw]*dat1$f
    dat1_tb <- data.table(dat1, key=t)
    names(dat1_tb)[ncol(dat1_tb)]="weight"
    names(dat1_tb)[names(dat1_tb)==t]="t_tb"
    sum_wt.t <- as.data.frame(dat1_tb[, sum(weight), by=list(t_tb)])
    #sum_wt.t = aggregate(dat1[,pw]*dat1$f, by=list(dat1[,t]), sum)
    sum_wt.t = sum_wt.t[match(unique(dat1[, t]), sum_wt.t[,1]),]
    unq_t.d_indx = d_indx[!duplicated(dat1[, t])]
  }else{
    sum_wt.t=cbind(dat1[, t], dat1$f*dat1[,pw])
    unq_t.d_indx = d_indx
  }
  
  #U_w_2.k.out =NULL
  if(is.null(post)){
    U_w_2 = 0 # Ai 
    for(j in 1:nrow(sum_wt.t)){
      k=unq_t.d_indx[j]
      U_w_2.k = sum_wt.t[j, 2]*(c((dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
                                  outer(c((dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
      U_w_2 = U_w_2.k+U_w_2
    }
    
    #for(k in d_indx){
    #  U_w_2.k = dat$f[k]*dat[k,pw]*(c((dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
    #                         outer(c((dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_2 = U_w_2.k+U_w_2
    #  U_w_2.k.out = cbind(U_w_2.k.out, U_w_2.k[,1])
    #  #print(k)
    #}
    Ui_pw= as.matrix(dat[,d]*(x.mtrx- H)-U_w_2)
  }else{
    if(nrow(f_w)==length(d_indx)){
      fw_indx=d_indx
    }else if(nrow(f_w)==nrow(dat)){
        fw_indx=1:nrow(dat)
    }else{
      stop("The partial derivative matrix of poststratification factor should be given for event cases or the whole sample!")
    }
    U_w_2 = 0 # Ai 
    U_w_4 = 0 # Bi
    #U_w_4.k.out =NULL
    for(j in 1:nrow(sum_wt.t)){
      k=unq_t.d_indx[j]
      U_w_2.k = sum_wt.t[j, 2]*(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
                                  outer(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
      U_w_2 = U_w_2.k+U_w_2
      #U_w_2.k.out = cbind(U_w_2.k.out, U_w_2.k[,1])
      U_w_4.k = sum_wt.t[j, 2]*(f_w_mtrx(f_w, (c(dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx)[fw_indx,])/dat$H_dnom[k]-
                                  outer(c(f_w_mtrx(f_w, c((dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)[fw_indx]))),
                                        as.numeric(as.matrix(dat$H_num)[k,]))/dat$H_dnom[k]/dat$H_dnom[k])
      U_w_4 = U_w_4.k+U_w_4
      #U_w_4.k.out = cbind(U_w_4.k.out,U_w_4.k[,1])
      #print(j)
    }
    #for(k in d_indx){
    #  U_w_2.k = dat$f[k]*dat[k,pw]*(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
    #                                  outer(c(dat$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_2 = U_w_2.k+U_w_2
    #  U_w_2.k.out=cbind(U_w_2.k.out, U_w_2.k[,1])
    #  U_w_4.k = dat$f[k]*dat[k,pw]*((f_w%*%as.matrix((c(dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx)[d_indx,]))/dat$H_dnom[k]-
    #                                  outer(c(f_w%*%c((dat[,pw]*(dat[,t]>=dat[,t][k])*dat$rel_hzd)[d_indx])),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    #  U_w_4 = U_w_4.k+U_w_4
    #  U_w_4.k.out = cbind(U_w_4.k.out,U_w_4.k[,1])
    #}
    if(length(fw_indx)== length(d_indx)){
      id_d = c(dat$id[d_indx], dat$id[dat[,d]==0])
      U_w_3 = rbind(f_w_mtrx(f_w, (dat[,pw]*(x.mtrx- H))[d_indx,]), matrix(0, n-sum(dat[,d]), p))[order(id_d),]
      U_w_4 = rbind(U_w_4, matrix(0, n-sum(dat[,d]), p))[order(id_d),]
    }else{
      U_w_3 = f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*(x.mtrx- H)))
      }
    Ui_pw= as.matrix(c(dat$f*dat[,d])*(x.mtrx- H)-U_w_2+U_w_3-U_w_4)
    
  }
  U_beta_1 = 0
  for(j in 1:nrow(sum_wt.t)){
    k=unq_t.d_indx[j]
    sum_j = t(as.matrix(x.mtrx*c(dat$pw_e*as.numeric(dat[,t]>=dat[,t][k]))))%*%as.matrix(x.mtrx)
    U_beta_1 =-sum_wt.t[j, 2]*sum_j/dat$H_dnom[k]+U_beta_1
  }
  #for(i in d_indx){
  #  sum_j = t(as.matrix(x.mtrx*c(dat$pw_e*as.numeric(dat[,t]>=dat[,t][i]))))%*%as.matrix(x.mtrx)
  #  U_beta_1 =-dat$f[i]*dat[i,pw]*sum_j/dat$H_dnom[i]+U_beta_1
  #}
  
  U_beta = U_beta_1 +t(c(dat[d_indx,pw])*H[d_indx,])%*%H[d_indx,]
  beta_pw = -Ui_pw%*%solve(U_beta)
  Delta_beta.pw = c(dat[,pw])*beta_pw
  
  var_beta_poisson.pw = t((1-pi.c_est)*Delta_beta.pw)%*%Delta_beta.pw
  var_beta_pps.pw = cov(Delta_beta.pw)*n
  U_list = list(Ui_pw  = Ui_pw,
                #U_w_2=U_w_2,U_w_3=U_w_3,U_w_4=U_w_4,
                U_beta = U_beta,
                x.mtrx = x.mtrx
  )
  #print("beta_pw.cox OK")
  return(list(U_list = U_list,
              beta_pw = beta_pw,
              Delta_beta.pw = Delta_beta.pw,
              var_beta_poisson.pw = var_beta_poisson.pw,
              var_beta_pps.pw = var_beta_pps.pw)
  )
}



lambda.ar_w = function(x.mtrx, dat, pw, rel_hzd, d="d", t="t", 
                       beta_wt=NULL, pw_list = NULL, post=NULL){
  n = nrow(dat)
  dat$id = 1:n
  dat$rel_hzd = rel_hzd
  Nt=NULL; Zt=NULL; Yt=NULL; t_unq = NULL
  #dat = dat[order(dat$t),]
  dat1 = dat[dat[,d]==1,]
  lambda_dat = data.frame(t = unique(dat1[,t]))
  names(lambda_dat)=t
  dat$Yi_t = outer(dat[,t], lambda_dat[,t], FUN=">=")
  colnames(dat$Yi_t)=paste("t=", lambda_dat[,t], sep="")
  row.names(dat$Yi_t) = row.names(dat)
  dat$Ii_t = outer(dat[,t], lambda_dat[,t], FUN="==")
  colnames(dat$Ii_t)=paste("t=", lambda_dat[,t], sep="")
  row.names(dat$Ii_t) = row.names(dat)

  for (i in 1:length(lambda_dat[,t])){
    Nt_i = sum((dat[,pw]*dat[,d])*dat$Ii_t[,i])
    Zt_i = sum((dat[,pw]*dat$rel_hzd)*dat$Yi_t[,i])
    Yt_i = sum(dat[,pw]*dat$Yi_t[,i])
    t_unq = c(t_unq, lambda_dat[,t][i])
    Nt = c(Nt, Nt_i)
    Zt = c(Zt, Zt_i)
    Yt = c(Yt, Yt_i)
  }  
  #for (t_i in lambda_dat$t){
  #  Nt_i = sum((dat[,pw]*dat$d)[dat$t==t_i])
  #  Zt_i = sum((dat[,pw]*dat$rel_hzd)[dat$t>=t_i])
  #  Yt_i = sum(dat[,pw]*(dat$t>=t_i))
  #  t_unq = c(t_unq, t_i)
  #  Nt = c(Nt, Nt_i)
  #  Zt = c(Zt, Zt_i)
  #  Yt = c(Yt, Yt_i)
  #}
  lambda_dat$Nt = Nt
  lambda_dat$Zt = Zt
  lambda_dat$Yt = Yt
  lambda_dat$lambda = Nt/Zt
  lambda_dat$ar = 1-Yt/Zt

  if(is.null(beta_wt)){
    return(list(lambda = lambda_dat$lambda,
                ar = lambda_dat$ar,
                u = lambda_dat[,t]))
    stop()
    
  }
  if(is.null(pw_list)){
    if(is.null(post)){
      Nt_w = c(dat[,d])*dat$Ii_t
      Zt_w = c(dat$rel_hzd)*dat$Yi_t+beta_wt%*%t(t(c(dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)
      Yt_w = dat$Yi_t
    }else{
      dat$f = post$f
      f_w = post$f_w
      dat[,pw] = dat[,pw]/dat$f
      d_indx = (dat[,d]==1)
      if(nrow(f_w)==sum(d_indx)){
        id_d = c(dat$id[d_indx], dat$id[1-d_indx])
        f_w_0 = matrix(0, n-sum(dat[,d]), length(t_unq))
        Nt_w = c(dat$f*dat[,d])*dat$Ii_t + rbind(f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*dat$Ii_t)[d_indx,]),
                                                 f_w_0)[order(id_d),]
        Zt_w = c(dat$f*dat$rel_hzd)*dat$Yi_t+beta_wt%*%t(t(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)+
          rbind(f_w_mtrx(f_w, (c(dat[,pw]*dat$rel_hzd)*dat$Yi_t)[d_indx,]), f_w_0)[order(id_d),]
        Yt_w = c(dat$f)*dat$Yi_t+rbind(f_w_mtrx(f_w, (c(dat[,pw])*dat$Yi_t)[d_indx,]), 
                                       f_w_0)[order(id_d),]
      }else if(nrow(f_w)==length(d_indx)){
        Nt_w = c(dat$f*dat[,d])*dat$Ii_t + f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*dat$Ii_t))
        Zt_w = c(dat$f*dat$rel_hzd)*dat$Yi_t+beta_wt%*%t(t(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)+
               f_w_mtrx(f_w, (c(dat[,pw]*dat$rel_hzd)*dat$Yi_t))
        Yt_w = c(dat$f)*dat$Yi_t+f_w_mtrx(f_w, (c(dat[,pw])*dat$Yi_t))
      }else{
        stop("The partial derivative matrix of poststratification factor should be given for event cases or the whole sample!")
      }
    }
  }
 if(!is.null(pw_list)){
   pw_gamma = pw_list$pw_gamma
   gamma_wt = pw_list$gamma_wt
   n_s = nrow(gamma_wt) - nrow(pw_gamma)
   if(is.null(post)){
     Nt_w = rbind(c(dat[,pw]*dat[,d])*dat$Ii_t, matrix(0, n_s, nrow(lambda_dat)))+
            gamma_wt%*%(t(pw_gamma)%*%(c(dat[,d])*dat$Ii_t))
     Zt_w = rbind(dat[,pw]*c(dat$rel_hzd)*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
            gamma_wt%*%(t(pw_gamma)%*%(c(dat$rel_hzd)*dat$Yi_t))+
            beta_wt%*%t(t(c(dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)
     Yt_w = rbind(dat[,pw]*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
            gamma_wt%*%(t(pw_gamma)%*%dat$Yi_t)
   }else{
     dat$f = post$f
     f_w = post$f_w
     dat[,pw] = dat[,pw]/dat$f
     d_indx = (dat[,d]==1)
     if(nrow(f_w)==sum(d_indx)){
       id_d = c(dat$id[d_indx], dat$id[1-d_indx])
       f_w_0 = matrix(0, n-sum(dat[,d]), length(t_unq))
       Nt_w = rbind(c(dat$f*dat[,pw]*dat[,d])*dat$Ii_t, matrix(0, n_s, nrow(lambda_dat)))+
         gamma_wt%*%(t(pw_gamma)%*%(c(dat$f*dat[,d])*dat$Ii_t+
                                      rbind(f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*dat$Ii_t)[d_indx,]), 
                                            f_w_0)[order(id_d),]))
       Zt_w = rbind(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
         gamma_wt%*%(t(pw_gamma)%*%(c(dat$f*dat$rel_hzd)*dat$Yi_t+
                                      rbind(f_w_mtrx(f_w, (c(dat[,pw]*dat$rel_hzd)*dat$Yi_t)[d_indx,]), 
                                            f_w_0)[order(id_d),]))+
         beta_wt%*%t(t(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)
       Yt_w = rbind(c(dat$f*dat[,pw])*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
         gamma_wt%*%(t(pw_gamma)%*%(c(dat$f)*dat$Yi_t+
                                      rbind(f_w_mtrx(f_w, (c(dat[,pw])*dat$Yi_t)[d_indx,]), 
                                            f_w_0)[order(id_d),]))
     }else if(nrow(f_w)==length(d_indx)){
       Nt_w = rbind(c(dat$f*dat[,pw]*dat[,d])*dat$Ii_t, matrix(0, n_s, nrow(lambda_dat)))+
              gamma_wt%*%(t(pw_gamma)%*%(c(dat$f*dat[,d])*dat$Ii_t+
                                         f_w_mtrx(f_w, (c(dat[,pw]*dat[,d])*dat$Ii_t))))
       Zt_w = rbind(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
              gamma_wt%*%(t(pw_gamma)%*%(c(dat$f*dat$rel_hzd)*dat$Yi_t+
                                         f_w_mtrx(f_w, (c(dat[,pw]*dat$rel_hzd)*dat$Yi_t))))+
              beta_wt%*%t(t(c(dat$f*dat[,pw]*dat$rel_hzd)*dat$Yi_t)%*%x.mtrx)
       Yt_w = rbind(c(dat$f*dat[,pw])*dat$Yi_t, matrix(0, n_s, nrow(lambda_dat)))+
              gamma_wt%*%(t(pw_gamma)%*%(c(dat$f)*dat$Yi_t+
                                         f_w_mtrx(f_w, (c(dat[,pw])*dat$Yi_t))))
     }else{
       stop("The partial derivative matrix of poststratification factor should be given for event cases or the whole sample!")
     }
   }

  }
    lambda_wt = t((lambda_dat$Zt)^-1*(t(Nt_w)-lambda_dat$lambda*t(Zt_w)))
    colnames(lambda_wt) = lambda_dat[,t]
    ar_wt = -t((lambda_dat$Zt)^-1*(t(Yt_w)-(1-lambda_dat$ar)*t(Zt_w)))
    colnames(ar_wt) = lambda_dat[,t]
    return(list(lambda = lambda_dat$lambda,
                ar = lambda_dat$ar,
                lambda_wt = lambda_wt,
                ar_wt = ar_wt,
                u = lambda_dat[,t])
           )
}
#t_star=c(1:5)
Lambda_w = function(lambda, lambda_wt=NULL, u, t_star){
  lambda    = lambda[order(u)]
  u = sort(u)
  t_diff = outer(u, t_star, FUN="-")
  t_diff[t_diff>0]=-10^6
  Lambda = cumsum(lambda)[apply(t_diff, 2, which.max)]
  if(is.null(lambda_wt)){
    return(list(Lambda = Lambda))
  }
  lambda_wt = lambda_wt[,order(u)]
  Lambda_wt = rowCumsums(lambda_wt)[,apply(t_diff, 2, which.max)]
  colnames(Lambda_wt) = paste("Lambda(t=", t_star, ")", sep="")
  return(list(Lambda = Lambda,
              Lambda_wt = Lambda_wt))
  
}


LambdaG_w = function(ar, ar_wt=NULL, u, t="t", lambda_star, t_star){
  ar.t = data.frame(ar = ar, t = u)
  ar.t = ar.t[order(ar.t[,t]),]

  lambda_star = lambda_star[lambda_star[,t]<=max(u),]
  bs_h_dat = merge(lambda_star, ar.t, all=T)
  final_t = nrow(bs_h_dat)
  cmp_t_indx = c(1:final_t)[!is.na(bs_h_dat$ar)]
  rep_time = cmp_t_indx-c(0, cmp_t_indx[-length(cmp_t_indx)])
  bs_h_dat$ar_cmp = rep(bs_h_dat$ar[cmp_t_indx], rep_time)
  
  t_diff = outer(bs_h_dat[,t], t_star, FUN="-")
  t_diff[t_diff>0]=-10^6
  LambdaG = cumsum(bs_h_dat$lambda_star*(1-bs_h_dat$ar_cmp))[apply(t_diff, 2, which.max)]
  if(is.null(ar_wt)){
    return(list(LambdaG = LambdaG))
    stop()
  }
    
  ar_wt = ar_wt[,order(u)]
  ar_wt.t = t(ar_wt)
  bs_h_dat$ar_wt.cmp = ar_wt.t[rep(1:nrow(ar_wt.t), rep_time),]
  
  LambdaG_wt = -t(colCumsums(bs_h_dat$lambda_star*bs_h_dat$ar_wt.cmp)[apply(t_diff, 2, which.max),])
  
  return(list(LambdaG = LambdaG, LambdaG_wt = LambdaG_wt))
}



absR_w = function(beta_est, Lambda, x0, beta_w=NULL, Lambda_w=NULL){
  n_beta = length(beta_est)
  x0 = matrix(x0, ncol=n_beta)
  beta_est = as.matrix(beta_est)
  rr = exp(x0%*%beta_est)
  absR = 1-exp(-outer(c(Lambda), c(rr)))
  if((!is.null(beta_w))&(!is.null(Lambda_w))){
    if((length(Lambda)>1)&(nrow(x0)==1)){
      n_t = length(Lambda)
      n = nrow(beta_w)
      beta_w_x0 = matrix(rep(beta_w%*%t(x0), n_t), n, n_t)
      absR_w = t(c((1-absR)*c(rr))*(Lambda*t(beta_w_x0)+t(Lambda_w)))
    }
    if((length(Lambda)==1)&(ncol(x0)>=1)){
      absR_w = t(c((1-absR)*rr)*t(Lambda*(beta_w%*%t(x0))+Lambda_w))
    }

    return(list(absR = absR,
                absR_w = absR_w))
  }else(return(list(absR = absR)))
}



post.f = function(samp, wt0, post.h, Nh, f_w=T, Large=F){
  n=length(samp[,wt0]);n
  if(length(post.h)==1){
    if(post.h%in%names(samp)){
      post.var = as.factor(samp[,post.h])
    }else{stop("Poststratification variable name must be included in the sample")}
  }else if(length(post.h)==n){
    post.var = as.factor(post.h)
  }else{stop("Poststratification variable must have the same length with the sample size")}
  post.var_num = as.numeric(post.var)
  id.h = model.matrix(~ post.var-1) # n * G groups
  Nh.hat = colSums(as.matrix(id.h*samp[,wt0]))  # sample totals in each group
  f = id.h%*%(Nh/Nh.hat); dim(f) # a vector of n of each Ng/Ng.hat
  if(f_w){
    if(!Large){
      f_w=-c(as.matrix(id.h)%*%as.matrix(Nh/Nh.hat/Nh.hat))*(id.h%*%t(id.h))
    }else{
      f_w=data.frame(f_w=-c(as.matrix(id.h)%*%as.matrix(Nh/Nh.hat/Nh.hat)), post.var = post.var_num)
      warning("R vector memory exhausted. Only f_w values are saved.")
    }
    
    return(list(f = c(f), f_w = f_w))    
  }else{
    return(list(f = c(f))) 
  }
}


#samp=data.frame(x1=1:5, x2=c(5,2,7,3,1), wt0=c(5:1))
#aux.tot = c("Intercept"=16, "x1"=60, "x2"=69)
#wt0="wt0"
greg.f = function(samp, wt0, aux.mtx=NULL, aux.tot, f_w=T){
  n=length(samp[,wt0]);n
  if(!is.null(aux.mtx)) samp = cbind(samp, aux.mtx)
  ds.wt = svydesign(ids=~1, data=samp, weights=as.formula(paste0("~",wt0)))
  calib.fm = as.formula(paste0("~", paste0(names(aux.tot)[-1],collapse="+")))
  V.hat = c(sum(samp[,wt0]), svytotal(calib.fm,ds.wt))
  V = aux.tot
  v.mtx = model.matrix(calib.fm,samp)
  vWv_inv = solve(t(samp[,wt0]*v.mtx)%*%v.mtx)
  f = c(t(1+(V-V.hat)%*%vWv_inv%*%t(v.mtx)))
  if(f_w){
    f_w1  = (V-V.hat)%*%vWv_inv
    vWv_w = lapply(1:n, function(i) outer(v.mtx[i,], v.mtx[i,]))
    f_w2  = vWv_inv%*%t(v.mtx)
    f_w = -sapply(1:n, function(i) f_w1%*%vWv_w[[i]]%*%f_w2)-
          (v.mtx)%*%vWv_inv%*%t(v.mtx)
    return(list(f = c(f), f_w = f_w))
  }else{return(list(f=c(f)))}
  
}





#######################################################################################
# cum_bsln_hzd is the FUNCTION to estimate cumulative baseline hazard                 #
# INPUT                                                                               #
#  dat     - a data frame of the sample that is used to estimate hazard               #
#  wt      - name of the sample weight variable in dat                                #
#  rel_hzd - a vector of relative hazard                                              #
#  t_star  - a vector of time point(s) which the baseline hazard is estimated for     #
# OUPUT                                                                               #
#  cum_hzd - a vector of estimates of cumulative baseline hazard at given time points #
#######################################################################################
cum_bsln_hzd = function(dat, wt, rel_hzd, t = "t", d = "d", t_star){
  dat$rel_hzd = rel_hzd
  dat = dat[order(dat[,t], -dat[,d]),]
  dat$w_e = dat[,wt]*dat$rel_hzd
  dat$dnom = rev(cumsum(rev(dat$w_e)))
  dat1 = dat[dat[,d]==1,c(t, d, wt, "dnom")]
  # Find ties in the dataset
  dup = duplicated(dat1[,t])
  if(sum(dup)>0){
    dnom_dat = dat1[!dup, c(t, "dnom")]
    dup_t = unique(dat1[,t][dup])
    dup_dat1 = dat1[dat1[,t]%in%dup_t, c(t, wt)]
    dup_de = aggregate(dup_dat1[,wt], by=list(dup_dat1[,t]), FUN=sum)
    names(dup_de) = c(t, wt)
    unq_dat1 = dat1[!dat1[,t]%in%dup_t, c(t, wt)]
    num_dat = rbind(unq_dat1, dup_de)
    num_dat = num_dat[order(num_dat[,t]),]
    #sum(dnom[,t]!=num_dat[,t])
    t_diff = outer(num_dat[,t], t_star, FUN="-")
    t_diff[t_diff>0]=-10^6
    lambda.u = num_dat[,wt]/dnom_dat$dnom
    cum_hzd = cumsum(lambda.u)[apply(t_diff, 2, which.max)]
  }else{
    t_diff = outer(dat1[,t], t_star, FUN="-")
    t_diff[t_diff>0]=-10^6
    lambda.u = dat1[,wt]/dat1$dnom
    cum_hzd = cumsum(lambda.u)[apply(t_diff, 2, which.max)]
  }
  return(list(cum_hzd  = cum_hzd,
              lambda.u = lambda.u, 
              u = unique(dat1[,t]))
  )
}
#######################################################################################
# gail_h is the FUNCTION to estimate cumulative baseline hazard using Gail's method   #
# INPUT                                                                               #
#  pop     - a data frame of the finite population                                    #
#  samp    - a data frame of the sample that is used to estimate hazard               #
#  wt      - name of the sample weight variable in samp                               #
#  rel_hzd - a vector of relative hazard                                              #
#  t_star  - a vector of time point(s) which the baseline hazard is estimated for     #
# OUPUT                                                                               #
#  cum_hzd - a vector of estimates of cumulative baseline hazard at given time points #
#######################################################################################
gail_h = function(pop=NULL, lambda_star=NULL, samp, wt, rel_hzd, t = "t", d = "d", t_star){
  if(is.null(lambda_star)){
    # lambda_star obtained from the finite population
    pop = pop[order(pop[,t], -pop[,d]),]
    pop$dnom = rev(cumsum(rev(pop$w)))
    pop1 = pop[pop[,d]==1,c(t, d, "w", "dnom")]
    ties_pop1 = duplicated(pop1[,t])
    if(sum(ties_pop1)>0){
      dnom_pop1 = pop1[!ties_pop1, c(t, "dnom")]
      dup_t = unique(pop1[,t][ties_pop1])
      dup_pop1 = pop1[pop1[,t]%in%dup_t, c(t, "w")]
      dedup_pop1 = aggregate(dup_pop1[,"w"], by=list(dup_pop1[,t]), FUN=sum)
      names(dedup_pop1) = c(t, "w")
      unq_pop1 = pop1[!pop1[,t]%in%dup_t, c(t, "w")]
      num_pop1 = rbind(unq_pop1, dedup_pop1)
      num_pop1 = num_pop1[order(num_pop1[,t]),]
      sum(dnom_pop1[,t]!=num_pop1[,t])
      lambda_star = data.frame(t = dnom_pop1[,t], 
                               lambda_star = c(num_pop1$w/dnom_pop1$dnom))
    }else{
      lambda_star = data.frame(t = pop1[,t], 
                               lambda_star = pop1$w/pop1$dnom)
    }
  }else{
  	lambda_star = data.frame(t = t_star, 
                             lambda_star = lambda_star)
  }

  # 1-AR
  samp$rel_hzd = rel_hzd
  samp = samp[order(samp[,t], -samp[,d]),]
  samp$w_e = samp[,wt]*samp$rel_hzd
  samp$dnom = rev(cumsum(rev(samp$w_e)))
  samp$num = rev(cumsum(rev(samp[,wt])))
  ties_samp = duplicated(samp[,t])
  dedup_samp = samp[(!ties_samp)&(samp[,d]==1),]
  samp_1_ar = data.frame(t = dedup_samp[,t], 
                         samp_1_ar = dedup_samp$num/dedup_samp$dnom)
  
  # cumulative attributable risk
  t_diff_ar = outer(samp_1_ar[,t], t_star, FUN="-")
  t_diff_ar[t_diff_ar>0]=-10^6
  est_1_ar = samp_1_ar$samp_1_ar[apply(t_diff_ar, 2, which.max)]
  bs_h_dat = merge(lambda_star, samp_1_ar, all=T)
  final_t = max(which(!is.na(bs_h_dat$samp_1_ar)))
  bs_h_dat = bs_h_dat[1:final_t,]
  
  cmp_t_indx = c(1:final_t)[!is.na(bs_h_dat$samp_1_ar)]
  rep_time = cmp_t_indx-c(0, cmp_t_indx[-length(cmp_t_indx)])
  
  bs_h_dat$samp_1_ar_cmp = rep(bs_h_dat$samp_1_ar[cmp_t_indx], rep_time)
  
  t_diff = outer(bs_h_dat[,t], t_star, FUN="-")
  t_diff[t_diff>0]=-10^6
  cum_hzd = cumsum(bs_h_dat$lambda_star*bs_h_dat$samp_1_ar_cmp)[apply(t_diff, 2, which.max)]
  return(list(lambda_star = lambda_star, samp_1_ar = est_1_ar, cum_hzd = cum_hzd))
}

U2 = function(beta_est, surv.fit, lambda_star, wt=NULL){
  var.names = all.vars(surv.fit$formula)
  t = var.names[1]
  d = var.names[2]
  dat = eval(surv.fit$call$data)
  x.mtrx = model.matrix(surv.fit)
  n = surv.fit$n
  if(is.null(dat)) dat = eval(surv.fit$survey.design$call$data)
  rel_hzd = exp(x.mtrx%*%beta_est)
  dat$rel_hzd = rel_hzd
  if(!is.null(surv.fit$survey.design)){
    dat$wt = 1/surv.fit$survey.design$prob
  }else if(!is.null(wt)){
    dat$wt=wt
  }else{
      dat$wt=1
  }

  dat$w_e = dat$wt*dat$rel_hzd
  dat$id = 1:n
  
  x.mtrx = as.matrix(x.mtrx[order(dat[,t]),])
  dat = dat[order(dat[,t]),]
  dat$H_num = as.matrix(cumsum(as.data.frame(c(dat$w_e)*x.mtrx)[n:1,]))[n:1,] #H1
  x.mtrx = as.matrix(x.mtrx[order(dat[,t], -dat[,d]),])
  dat = dat[order(dat[,t], -dat[,d]),]
  dat$dnom = rev(cumsum(rev(dat$w_e)))
  dat$num = rev(cumsum(rev(dat$wt)))
  ties_dat = duplicated(dat[,t])
  dedup_dat = dat[(!ties_dat)&(dat[,d]==1),]
  dat_1_ar = data.frame(t = dedup_dat[,t], 
                        dat_1_ar = dedup_dat$num/dedup_dat$dnom,
                        wt = dedup_dat$wt
  )
  names(lambda_star)[1]=t
  bs_h_dat = merge(lambda_star, dat_1_ar, all.y=T)
  U1 = colSums(matrix((dat$wt*x.mtrx)[dat[,d]==1,], ncol=length(beta_est)))
  U2 = colSums(matrix((bs_h_dat$lambda_star*dat_1_ar$wt*bs_h_dat$dat_1_ar)*dedup_dat$H_num, ncol=length(beta_est)))
  U = U1-U2
  sum(abs(U))
}

U = function(beta_est, surv.fit){
  var.names = all.vars(surv.fit$formula)
  t = var.names[1]
  d = var.names[2]
  dat = eval(surv.fit$call$data)
  x.mtrx = model.matrix(surv.fit)
  n = surv.fit$n
  if(is.null(dat)) dat = eval(surv.fit$survey.design$call$data)
  rel_hzd = exp(x.mtrx%*%beta_est)
  dat$rel_hzd = rel_hzd
  dat$wt=1;
  if(!is.null(surv.fit$weights)) dat$wt = surv.fit$weights*surv.fit$n
  dat$w_e = dat$wt*dat$rel_hzd
  dat$id = 1:n
  x.mtrx = as.matrix(x.mtrx[order(dat[,t]),])
  dat = dat[order(dat[,t]),]
  dat$H_dnom = rev(cumsum(rev(dat$w_e)))# H2
  dat$H_num = as.matrix(cumsum(as.data.frame(c(dat$w_e)*x.mtrx)[n:1,]))[n:1,] #H1
  ties = duplicated(dat[,t])
  if(sum(ties)>0){
    H_uniq = dat[!ties,c(t, "H_dnom", "H_num")]
    H_dat = H_uniq[rep(1:nrow(H_uniq), table(dat[,t])),]
    #h_dat[,t]==dat[,t]
    dat[, c("H_dnom", "H_num")] = H_dat[, c("H_dnom", "H_num")]
  }
  x.mtrx = as.matrix(x.mtrx[order(dat$id), ])
  dat = dat[order(dat$id), ]
  
  H = as.matrix(dat$H_num/dat$H_dnom)
  U = colSums(dat[,d]*(x.mtrx- H))
  sum(abs(U))
}

lambda_star.pop = function(pop, t="t", d="d", wt=NULL, t_star=NULL){
  if(is.null(wt)) {
    wt="wt"
    pop[,wt]=1}
  pop = pop[order(pop[,t], -pop[,d]),]
  pop$dnom = rev(cumsum(rev(pop[,wt])))
  pop1 = pop[pop[,d]==1,c(t, d, wt, "dnom")]
  ties_pop1 = duplicated(pop1[,t])
  if(sum(ties_pop1)>0){
    dnom_pop1 = pop1[!ties_pop1, c(t, "dnom")]
    dup_t = unique(pop1[,t][ties_pop1])
    dup_pop1 = pop1[pop1[,t]%in%dup_t, c(t, wt)]
    dedup_pop1 = aggregate(dup_pop1[,wt], by=list(dup_pop1[,t]), FUN=sum)
    names(dedup_pop1) = c(t, wt)
    unq_pop1 = pop1[!pop1[,t]%in%dup_t, c(t, wt)]
    num_pop1 = rbind(unq_pop1, dedup_pop1)
    num_pop1 = num_pop1[order(num_pop1[,t]),]
    sum(dnom_pop1[,t]!=num_pop1[,t])
    lambda_star = data.frame(t = dnom_pop1[,t], 
                             lambda_star = c(num_pop1[,wt]/dnom_pop1$dnom),
                             event_set = num_pop1[,wt],
                             risk_set = dnom_pop1$dnom)
  }else{
    lambda_star = data.frame(t = pop1[,t], 
                             lambda_star = pop1[,wt]/pop1$dnom,
                             event_set = pop1[,wt],
                             risk_set = pop1$dnom)
  }
  lambda_star = lambda_star[order(lambda_star$t),]
  if(!is.null(t_star)){
    t_star=sort(t_star)
    
    lambda_star$t_int = cut(lambda_star$t, breaks=c(0,t_star, max(lambda_star$t)+1), include.lowest = T)
    lambda_tb <- data.table(lambda_star, key="t_int")
    event_set <- as.data.frame(lambda_tb[, sum(event_set), by=list(t_int)])
    keep = !rev(duplicated(rev(lambda_star$t_int)))
    lambda_star_sub = lambda_star[keep,]
    lambda_star_sub$event_set=event_set[,2]
    lambda_star_sub$lambda_star=lambda_star_sub$event_set/lambda_star_sub$risk_set
    lambda_star = lambda_star_sub[,c(1:4)]
  }
  lambda_star
}

var_Delta = function(Delta, R=NULL){
  Delta.dat = as.data.frame(Delta)
  names(Delta.dat) = paste0("Delta.", 1:ncol(Delta))
  ds = svydesign(ids=~1,strata = R, weight=1, data = Delta.dat)
  var.pps = diag(vcov(svytotal(as.formula(paste0("~", paste0(names(Delta.dat), collapse="+"))), ds)))
  var.pps
}

