est_out = function(cox.fit, samp, pw.name, d="d", t="t", betas=NULL, x0=NULL){
  if(is.null(betas)) betas=cox.fit$coeff
  x.mtrx = model.matrix(cox.fit)
  rel_hzd = c(exp(x.mtrx%*%betas))
  lmda = lambda.ar_w(x.mtrx = x.mtrx, rel_hzd = rel_hzd, d=d, t=t,
                     dat = samp, pw = pw.name)  
  Lmda  =  Lambda_w(lambda = lmda$lambda, u = lmda$u,t_star = t_star)
  LmdaG =  LambdaG_w(ar=lmda$ar, u=lmda$u,
                     #lambda_star = gail_pop_out$lambda_star, 
                     lambda_star = lambda_star, 
                     t_star=t_star)
  absR  = absR_w(beta_est = betas, Lambda   = Lmda$Lambda,   x0 =x0)$absR
  absRG = absR_w(beta_est = betas, Lambda   = LmdaG$LambdaG, x0 =x0)$absR
  #absR  = absR_w(beta_est = betas, Lambda   = Lmda$Lambda,   x0 =x0[,3])
  #absRG = absR_w(beta_est = betas, Lambda   = LmdaG$LambdaG, x0 =x0[,3])
  
  return(list(Lambda      = Lmda$Lambda,   
              LambdaG     = LmdaG$LambdaG, 
              absR        = absR,     
              absRG       = absRG    
  ))
}

calib_est = function(var.t, var.d, var.wt, com_dat, a_cv=NULL){
  surv.fm = paste0("Surv(",var.t,",",var.d, ") ~ x1+x2+x3")
  samp.c=com_dat[com_dat$trt==1,]
  if(!var.wt=="wt"){
    samp.c[,var.d] = samp.c[,var.d]/a_cv[1]
  } else{
    samp.c$wt = samp.c$wt*2
  }
  samp.s=com_dat[com_dat$trt==0,]
  samp.s$wt = samp.s$wt*2
  ds.comb = svydesign(ids=~1, strata=~trt, data=com_dat, 
                      weights=as.formula(paste0("~",var.wt)))
  aux.fit = svycoxph(as.formula(surv.fm), design = ds.comb)
  rr = c(exp(model.matrix(aux.fit)%*%aux.fit$coefficients))
  beta_w.aux = beta_pw.cox(x.mtrx = model.matrix(aux.fit), 
                           rel_hzd = rr, 
                           dat = com_dat, pw = var.wt, t=var.t, d=var.d)
  aux.mtx = data.frame(cnt = 1, Delta_beta = beta_w.aux$beta_pw[com_dat$trt==1,],
                       t.rr = samp.c[,var.t]*(rr[com_dat$trt==1]))
  #names(aux.mtx)[1]="(Intercept)"
  aux.tot = c(N, N1, colSums(com_dat[,var.wt]*cbind(matrix(0,nrow(com_dat),n_beta), 
                                                    com_dat[,var.t]*rr)))
  names(aux.tot) = c("(Intercept)", "d", paste0("Delta_beta.x", c(1:n_beta)),"t.rr")
  ds.wt = svydesign(ids=~1, data=cbind(samp.c, aux.mtx), weights=as.formula(paste0("~",var.wt)))
  ##############pop tot and d prev
  #calib.rr = weights(calibrate(ds.wt, as.formula(paste0("~", paste0(names(aux.tot)[-c(1,3,length(aux.tot))],collapse="+"))), 
  #                               aux.tot[-c(3,length(aux.tot))], calfun="linear"))
  greg_beta = greg.f(samp=samp.c, wt0=var.wt, aux.mtx=aux.mtx, aux.tot=aux.tot[-length(aux.tot)], f_w=F)
  #greg_beta = greg.f(samp=samp.c, wt0=var.wt, aux.mtx=aux.mtx, aux.tot=aux.tot[c(1,2)], f_w=T)
  calib.rr = samp.c[,var.wt]*greg_beta$f
  calib.rr[calib.rr<0]=1e-5
  samp.c$calibwt.rr = calib.rr
  # Beta estimates 
  cox.c_calib = svycoxph(fm_fit.cox, design = svydesign(ids=~1, data=samp.c, weights=calib.rr))#;cox_out.1
  
  #calib.lamd = weights(calibrate(ds.wt, as.formula(paste0("~", paste0(names(aux.tot)[-c(1,3,1:n_beta+3)],collapse="+"))), 
  #                               aux.tot[c(1:2, length(aux.tot))], calfun="linear"))
  gw_lmd = greg.f(samp=samp.c, wt0=var.wt, aux.mtx=aux.mtx, aux.tot=aux.tot[c(1,2, length(aux.tot))], f_w=F)
  #gw_lmd = greg.f(samp=samp.c, wt0=var.wt, aux.mtx=aux.mtx, aux.tot=aux.tot[c(1,2)], f_w=T)
  calib.lamd=samp.c[,var.wt]*gw_lmd$f
  calib.lamd[calib.lamd<0]=1e-5
  
  est.calib   = est_out(cox.fit = cox.c_calib, pw.name = "calib.wt",
                        samp = data.frame(samp.c, calib.wt=calib.lamd))
  
  
  return(list(betas     = cox.c_calib$coefficients,
              ps.wts    = cbind(calib.rr, calib.lamd),
              Lambda    = est.calib$Lambda,  
              LambdaG   = est.calib$LambdaG, 
              absR      = est.calib$absR,    
              absRG     = est.calib$absRG   
  ))
  
}
jk_est_noninf = function(samp.c, samp.s, recal.wt=T){
  samp.c<<-samp.c
  n.mth=length(mth.names)
  n_s=nrow(samp.s)
  n_c=nrow(samp.c)
  beta_est    = matrix(0, n_beta, n.mth)
  Lambda_est  = matrix(0, length(t_star), n.mth)
  LambdaG_est = matrix(0, length(t_star), n.mth)
  absR_est    = array (0, c(length(t_star), n.mth, nrow(x0)))
  absRG_est   = array (0, c(length(t_star), n.mth, nrow(x0)))
  # Naive
  cox.c1 = coxph(fm_fit.cox, data = samp.c, robust=T, ties="breslow")
  beta_est[,1] = cox.c1$coefficients
  est = est_out(cox.fit=cox.c1, samp=samp.c, pw.name="w")
  Lambda_est[,1 ] = est$Lambda;  LambdaG_est[,1 ] = est$LambdaG
  absR_est  [,1,] = est$absR;    absRG_est  [,1,] = est$absRG
  # Cht
  cox.c = svycoxph(fm_fit.cox, design = svydesign(ids=~1, data=samp.c, weights=~wt))
  rel_hzd = exp(model.matrix(cox.c)%*%c(cox.c$coefficients))
  beta_est[,2] = cox.c$coefficients
  est = est_out(cox.fit=cox.c, samp = samp.c, pw.name="wt")
  Lambda_est[,2 ] = est$Lambda; LambdaG_est[,2 ] = est$LambdaG
  absR_est  [,2,] = est$absR;   absRG_est  [,2,] = est$absRG
  # Svy
  cox.s = svycoxph(fm_fit.cox, design = svydesign(ids=~1, data=samp.s, weights=~wt))
  beta_est[,3] = cox.s$coefficients
  est = est_out(cox.s, samp.s, "wt")
  Lambda_est[,3 ] = est$Lambda; LambdaG_est[,3 ] = est$LambdaG
  absR_est  [,3,] = est$absR;   absRG_est  [,3,] = est$absRG
  ## Combine samples
  samp.c$trt=1
  samp.s$trt=0
  com_dat = rbind(samp.c, samp.s)
  com_dat$wt=com_dat$wt/2
  for(i in 1:3){
    com_dat[,paste0("t_fit.",i)] = com_dat$t
    com_dat[,paste0("d_fit.",i)] = com_dat$d
    com_dat[com_dat$trt==0,paste0("t_fit.",i)] = samp.s[,paste0("t.imp.", i)]
    com_dat[com_dat$trt==0,paste0("d_fit.",i)] = samp.s[,paste0("d_tilde.", i)]
  }
  if(recal.wt){
    # IPSW.S
    a=n_s/N
    com_dat$w.s=com_dat$w
    com_dat$w.s[com_dat$trt==0]=samp.s$wt*a
    ps.ds = svydesign(ids =~1, strata=~trt, weight = ~w.s, data = com_dat)
    ps.model.par = ps.model.fit(ps.ds, fm_fit.s=fm_ps, a=n_s/N)
    samp.c$ipsw = 1/ps.model.par$model.par$pi.c_est
  }
  
  ds.c = svydesign(ids=~1, data=samp.c, weights=~ipsw, nest=TRUE)
  cox_out = svycoxph(fm_fit.cox, design = ds.c)
  beta_est[,4] = cox_out$coefficients
  est = est_out(cox.fit=cox_out, pw.name="ipsw", samp=samp.c)
  Lambda_est[,4 ] = est$Lambda; LambdaG_est[,4 ] = est$LambdaG
  absR_est  [,4,] = est$absR;   absRG_est  [,4,] = est$absRG
  # Calibration
  cv = (c(var(samp.c$ipsw)/mean(samp.c$ipsw^2), var(samp.s$wt)/mean(samp.s$wt^2))+1)/c(n_c, n_s)
  a_cv = (sum(samp.c$ipsw)+sum(samp.s$wt))*(1-cv/sum(cv))/c(sum(samp.c$ipsw), sum(samp.s$wt))/2
  #a_cv=c(.5, .5)
  com_dat$ipsw = c(samp.c$ipsw*a_cv[1], samp.s$wt*a_cv[2])
  for(i in 1:3){
    # Calibration using the true weights
    est = calib_est(var.t=paste0("t_tilde.",i), var.d=paste0("d_tilde.",i), var.wt="wt", com_dat=com_dat)  
    beta_est   [,5+(i-1)*4 ] = est$betas;    
    Lambda_est [,5+(i-1)*4 ] = est$Lambda;   LambdaG_est[,5+(i-1)*4 ] = est$LambdaG
    absR_est   [,5+(i-1)*4,] = est$absR;     absRG_est  [,5+(i-1)*4,] = est$absRG
    # Calibration+imputation
    #date()
    est = calib_est(var.t=paste0("t_fit.",i),var.d=paste0("d_fit.",i), var.wt="wt",com_dat=com_dat)  
    #date()
    beta_est   [,7+(i-1)*4 ] = est$betas; 
    Lambda_est [,7+(i-1)*4 ] = est$Lambda;   LambdaG_est[,7+(i-1)*4 ] = est$LambdaG
    absR_est   [,7+(i-1)*4,] = est$absR;     absRG_est  [,7+(i-1)*4,] = est$absRG
    
    ##########################################
    #           Method 1 Cmb-Calib           #
    ##########################################
    est = calib_est(var.t=paste0("t_tilde.",i),var.d=paste0("d_tilde.",i),var.wt="ipsw",com_dat=com_dat,a_cv=a_cv)  
    beta_est   [,6+(i-1)*4 ] = est$betas;
    Lambda_est [,6+(i-1)*4 ] = est$Lambda;   LambdaG_est[,6+(i-1)*4 ] = est$LambdaG
    absR_est   [,6+(i-1)*4,] = est$absR;     absRG_est  [,6+(i-1)*4,] = est$absRG
    est = calib_est(var.t=paste0("t_fit.",i),var.d=paste0("d_fit.",i),var.wt="ipsw",com_dat=com_dat,a_cv=a_cv)  
    beta_est   [,8+(i-1)*4 ] = est$betas;
    Lambda_est [,8+(i-1)*4 ] = est$Lambda;   LambdaG_est[,8+(i-1)*4 ] = est$LambdaG
    absR_est   [,8+(i-1)*4,] = est$absR;     absRG_est  [,8+(i-1)*4,] = est$absRG
    #print(i)
  }
  return(list(beta_est=beta_est, Lambda_est=Lambda_est, LambdaG_est=LambdaG_est,
              absR_est=absR_est, absRG_est=absRG_est))
  
}

jk_est_inf = function(samp.c, samp.s, recal.wt=T){
  samp.c<<-samp.c
  n.mth=length(mth.names)
  n_s=nrow(samp.s)
  n_c=nrow(samp.c)
  beta_est    = matrix(0, n_beta, n.mth)
  Lambda_est  = matrix(0, length(t_star), n.mth)
  LambdaG_est = matrix(0, length(t_star), n.mth)
  absR_est    = array (0, c(length(t_star), n.mth, nrow(x0)))
  absRG_est   = array (0, c(length(t_star), n.mth, nrow(x0)))
  # Naive
  cox.c1 = coxph(fm_fit.cox, data = samp.c, robust=T, ties="breslow")
  beta_est[,1] = cox.c1$coefficients
  est = est_out(cox.fit=cox.c1, samp=samp.c, pw.name="w")
  Lambda_est[,1]  = est$Lambda;  LambdaG_est[,1]  = est$LambdaG
  absR_est  [,1,] = est$absR;    absRG_est  [,1,] = est$absRG
  # Cht
  cox.c = svycoxph(fm_fit.cox, design = svydesign(ids=~1, data=samp.c, weights=~wt))
  rel_hzd = exp(model.matrix(cox.c)%*%c(cox.c$coefficients))
  beta_est[,2] = cox.c$coefficients
  est = est_out(cox.fit=cox.c, samp = samp.c, pw.name="wt")
  Lambda_est[,2 ] = est$Lambda; LambdaG_est[,2 ] = est$LambdaG
  absR_est  [,2,] = est$absR;   absRG_est  [,2,] = est$absRG
  # Svy
  cox.s = svycoxph(fm_fit.cox, design = svydesign(ids=~1, data=samp.s, weights=~wt))
  beta_est[,3] = cox.s$coefficients
  est = est_out(cox.s, samp.s, "wt")
  Lambda_est[,3 ] = est$Lambda; LambdaG_est[,3 ] = est$LambdaG
  absR_est  [,3,] = est$absR;   absRG_est  [,3,] = est$absRG
  ## Combine samples
  samp.c$trt=1
  samp.s$trt=0
  com_dat = rbind(samp.c, samp.s)
  com_dat$wt=com_dat$wt/2
  for(i in 1:3){
    com_dat[,paste0("t_fit.",i)] = com_dat$t
    com_dat[,paste0("d_fit.",i)] = com_dat$d
    com_dat[com_dat$trt==0,paste0("t_fit.",i)] = samp.s[,paste0("t.imp.", i)]
    com_dat[com_dat$trt==0,paste0("d_fit.",i)] = samp.s[,paste0("d_tilde.", i)]
  }
  
  # IPSW.S
  a=n_s/N
  com_dat$w.s=com_dat$w
  com_dat$w.s[com_dat$trt==0]=samp.s$wt*a
  ps.ds = svydesign(ids =~1, strata=~trt, weight = ~w.s, data = com_dat)
  a_cv=NULL
  for(i in 1:4){
    ipsw.name=paste0("ipsw.",i-1)
    ps.model.par = ps.model.fit(ps.ds, fm_fit.s=as.formula(fm_ps[i]), a=n_s/N)
    ipsw = 1/ps.model.par$model.par$pi.c_est
    ds.c = svydesign(ids=~1, data=samp.c, weights=ipsw, nest=TRUE)
    cox_out = svycoxph(fm_fit.cox, design = ds.c)
    beta_est[,3+i] = cox_out$coefficients
    est = est_out(cox.fit = cox_out, pw.name = "ipsw",  
                  samp = cbind(samp.c, ipsw=ipsw))
    Lambda_est[,3+i ] = est$Lambda;     LambdaG_est[,3+i ] = est$LambdaG
    absR_est  [,3+i,] = est$absR;       absRG_est  [,3+i,] = est$absRG
    if(recal.wt) samp.c[,ipsw.name]=ipsw
    
    cv = (c(var(samp.c[,ipsw.name])/mean(samp.c[,ipsw.name]^2), var(samp.s$wt)/mean(samp.s$wt^2))+1)/c(n_c, n_s)
    a_cv=rbind(a_cv, 
               (sum(samp.c[,ipsw.name])+sum(samp.s$wt))*(1-cv/sum(cv))/c(sum(samp.c[,ipsw.name]), 
                                                                         sum(samp.s$wt))/2)
    #a_cv=c(.5, .5)
    com_dat[,ipsw.name] = c(samp.c[,ipsw.name]*a_cv[i,1], samp.s$wt*a_cv[i,2])
  }
  # Calibration
  for(i in 1:3){
    # Calibration using the true weights
    est = calib_est(var.t=paste0("t_tilde.",i), var.d=paste0("d_tilde.",i), var.wt="wt", com_dat=com_dat)  
    beta_est   [,8+(i-1)*6 ] = est$betas
    Lambda_est [,8+(i-1)*6 ] = est$Lambda;   LambdaG_est[,8+(i-1)*6 ] = est$LambdaG
    absR_est   [,8+(i-1)*6,] = est$absR;     absRG_est  [,8+(i-1)*6,] = est$absRG
    # Calibration+imputation
    est = calib_est(var.t=paste0("t_fit.",i),var.d=paste0("d_fit.",i), var.wt="wt",com_dat=com_dat)  
    beta_est   [,11+(i-1)*6 ] = est$betas
    Lambda_est [,11+(i-1)*6 ] = est$Lambda;   LambdaG_est[,11+(i-1)*6 ] = est$LambdaG
    absR_est   [,11+(i-1)*6,] = est$absR;     absRG_est  [,11+(i-1)*6,] = est$absRG
    
    ##########################################
    #           Method 1 Cmb-Calib           #
    ##########################################
    est = calib_est(var.t=paste0("t_tilde.",i),var.d=paste0("d_tilde.",i),var.wt="ipsw.0",com_dat=com_dat,a_cv=a_cv[1,])  
    beta_est   [,9+(i-1)*6 ] = est$betas
    Lambda_est [,9+(i-1)*6 ] = est$Lambda;   LambdaG_est[,9+(i-1)*6 ] = est$LambdaG
    absR_est   [,9+(i-1)*6,] = est$absR;     absRG_est  [,9+(i-1)*6,] = est$absRG
    est = calib_est(var.t=paste0("t_tilde.",i),var.d=paste0("d_tilde.",i),var.wt=paste0("ipsw.",i),com_dat=com_dat,a_cv=a_cv[i+1,])  
    beta_est   [,10+(i-1)*6 ] = est$betas
    Lambda_est [,10+(i-1)*6 ] = est$Lambda;   LambdaG_est[,10+(i-1)*6 ] = est$LambdaG
    absR_est   [,10+(i-1)*6,] = est$absR;     absRG_est  [,10+(i-1)*6,] = est$absRG
    est = calib_est(var.t=paste0("t_fit.",i),var.d=paste0("d_fit.",i),var.wt="ipsw.0",com_dat=com_dat, a_cv=a_cv[1,])  
    beta_est   [,12+(i-1)*6 ] = est$betas
    Lambda_est [,12+(i-1)*6 ] = est$Lambda;   LambdaG_est[,12+(i-1)*6 ] = est$LambdaG
    absR_est   [,12+(i-1)*6,] = est$absR;     absRG_est  [,12+(i-1)*6,] = est$absRG
    est = calib_est(var.t=paste0("t_fit.",i),var.d=paste0("d_fit.",i),var.wt=paste0("ipsw.",i),com_dat=com_dat,a_cv=a_cv[i+1,])  
    beta_est   [,13+(i-1)*6 ] = est$betas
    Lambda_est [,13+(i-1)*6 ] = est$Lambda;   LambdaG_est[,13+(i-1)*6 ] = est$LambdaG
    absR_est   [,13+(i-1)*6,] = est$absR;     absRG_est  [,13+(i-1)*6,] = est$absRG
  }
  return(list(beta_est=beta_est, Lambda_est=Lambda_est, LambdaG_est=LambdaG_est,
              absR_est=absR_est, absRG_est=absRG_est))
  
}



jk_fun=function(samp.c, samp.s, m_jk, n_jk, recal.wt=T, sampling){
  n_c = nrow(samp.c)
  n_s = nrow(samp.s)
  tmp = runif(n_c)
  samp.c$jk_group = as.numeric(cut(tmp, quantile(tmp, prob=seq(0, 1, 1/m_jk)), include.lowest=T))
  #table(jk_group)
  tmp = runif(n_s)
  samp.s$jk_group = as.numeric(cut(tmp, quantile(tmp, prob=seq(0, 1, 1/n_jk)), include.lowest=T))
  #table(jk_group)
  beta_est=NULL;  Lambda_est=NULL;  LambdaG_est=NULL; 
  l.absR_est=NULL ; m.absR_est=NULL ; h.absR_est=NULL ; 
  l.absRG_est=NULL; m.absRG_est=NULL; h.absRG_est=NULL; 
  if(!recal.wt){
    samp.c$trt=1
    samp.s$trt=0
    com_dat = rbind(samp.c, samp.s)
    # IPSW.S
    a=n_s/N
    com_dat$w.s=com_dat$w
    com_dat$w.s[com_dat$trt==0]=samp.s$wt*a
    ps.ds = svydesign(ids =~1, strata=~trt, weight = ~w.s, data = com_dat)
    if(sampling=="noninf"){
      ps.model.par = ps.model.fit(ps.ds, fm_fit.s=fm_ps, a=n_s/N)
      samp.c$ipsw = 1/ps.model.par$model.par$pi.c_est
      samp.s$ipsw=samp.s$wt
    }
    if(sampling=="inf"){
      for(i in 1:4){
        ipsw.name=paste0("ipsw.",i-1)
        ps.model.par = ps.model.fit(ps.ds, fm_fit.s=as.formula(fm_ps[i]), a=n_s/N)
        samp.c[,ipsw.name] = 1/ps.model.par$model.par$pi.c_est
        samp.s[,ipsw.name]=samp.s$wt
      }
    }
  }
  for (k in 1:m_jk){
    # remove one psu at each replicate
    samp.c.k = samp.c[samp.c$jk_group!=k,]
    samp.c.k$w = m_jk/(m_jk-1)
    samp.c.k$wt= samp.c.k$wt*m_jk/(m_jk-1)
    if(sampling=="noninf"){
      if(!recal.wt){
        samp.c.k$ipsw = samp.c.k$ipsw*m_jk/(m_jk-1)
        jk.c.k = jk_est_noninf(samp.c=samp.c.k, samp.s=samp.s, recal.wt=F)        
      } else{
        jk.c.k = jk_est_noninf(samp.c=samp.c.k, samp.s=samp.s, recal.wt=T)        
      }
    } 
    if(sampling=="inf"){
      if(!recal.wt){
        for(i in 1:4){
          samp.c.k[,paste0("ipsw.",i-1)] = samp.c.k[,paste0("ipsw.",i-1)]*m_jk/(m_jk-1)
        }
        jk.c.k = jk_est_inf(samp.c=samp.c.k, samp.s=samp.s, recal.wt=F)        
      }else{
        jk.c.k = jk_est_inf(samp.c=samp.c.k, samp.s=samp.s, recal.wt=T)        
      }
    }    
    n.mth   = ncol(jk.c.k$beta_est)
    n_beta  = nrow(jk.c.k$beta_est)
    n_tstar = nrow(jk.c.k$Lambda_est)
    beta_est    = rbind(beta_est,    c(jk.c.k$beta_est))
    Lambda_est  = rbind(Lambda_est,  c(jk.c.k$Lambda_est))
    LambdaG_est = rbind(LambdaG_est, c(jk.c.k$LambdaG_est))
    l.absR_est    = rbind(l.absR_est ,   c(jk.c.k$absR_est [,,1]))
    l.absRG_est   = rbind(l.absRG_est,   c(jk.c.k$absRG_est[,,1]))
    m.absR_est    = rbind(m.absR_est ,   c(jk.c.k$absR_est [,,2]))
    m.absRG_est   = rbind(m.absRG_est,   c(jk.c.k$absRG_est[,,2]))
    h.absR_est    = rbind(h.absR_est ,   c(jk.c.k$absR_est [,,3]))
    h.absRG_est   = rbind(h.absRG_est,   c(jk.c.k$absRG_est[,,3]))
    print(paste0("c.",k))
  }
  for (k in 1:n_jk){
    # remove one psu at each replicate
    samp.s.k = samp.s[samp.s$jk_group!=k,]
    samp.s.k$w= n_jk/(n_jk-1)
    samp.s.k$wt= samp.s.k$wt*n_jk/(n_jk-1)
    if(!recal.wt){
      if(sampling=="noninf") jk.s.k = jk_est_noninf(samp.c=samp.c, samp.s=samp.s.k, recal.wt=F)
      if(sampling=="inf")    jk.s.k = jk_est_inf   (samp.c=samp.c, samp.s=samp.s.k, recal.wt=F)
    }else{
      if(sampling=="noninf") jk.s.k = jk_est_noninf(samp.c=samp.c, samp.s=samp.s.k, recal.wt=T)
      if(sampling=="inf")    jk.s.k = jk_est_inf   (samp.c=samp.c, samp.s=samp.s.k, recal.wt=T)
      
    }    
    beta_est    = rbind(beta_est,    c(jk.s.k$beta_est))
    Lambda_est  = rbind(Lambda_est,  c(jk.s.k$Lambda_est))
    LambdaG_est = rbind(LambdaG_est, c(jk.s.k$LambdaG_est))
    l.absR_est    = rbind(l.absR_est ,   c(jk.c.k$absR_est [,,1]))
    l.absRG_est   = rbind(l.absRG_est,   c(jk.c.k$absRG_est[,,1]))
    m.absR_est    = rbind(m.absR_est ,   c(jk.c.k$absR_est [,,2]))
    m.absRG_est   = rbind(m.absRG_est,   c(jk.c.k$absRG_est[,,2]))
    h.absR_est    = rbind(h.absR_est ,   c(jk.c.k$absR_est [,,3]))
    h.absRG_est   = rbind(h.absRG_est,   c(jk.c.k$absRG_est[,,3]))
    print(paste0("s.",k))
  }
  return(list(
    beta_est    = array(beta_est  ,  c(m_jk+n_jk, n_beta,  n.mth)),
    Lambda_est  = array(Lambda_est,  c(m_jk+n_jk, n_tstar, n.mth)),
    LambdaG_est = array(LambdaG_est, c(m_jk+n_jk, n_tstar, n.mth)),
    l.absR_est    = array(l.absR_est  ,  c(m_jk+n_jk, n_tstar, n.mth)),
    l.absRG_est   = array(l.absRG_est ,  c(m_jk+n_jk, n_tstar, n.mth)),
    m.absR_est    = array(m.absR_est  ,  c(m_jk+n_jk, n_tstar, n.mth)),
    m.absRG_est   = array(m.absRG_est ,  c(m_jk+n_jk, n_tstar, n.mth)),
    h.absR_est    = array(h.absR_est  ,  c(m_jk+n_jk, n_tstar, n.mth)),
    h.absRG_est   = array(h.absRG_est ,  c(m_jk+n_jk, n_tstar, n.mth))
  ))
  
}
