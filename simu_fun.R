########################################################################################
## simu_fun is a function to generate IPSW, PSAS, and KW weights                      ##
## INPUT                                                                              ##
##   chtsamp:   cohort sample data                                                    ##
##   svysamp:   survey sample data                                                    ##
##   svy_wt:    variable name of survey weights (should be included in "svysamp")     ##
##   y:         name of variable of interest (in cohort)                              ##
##   Formula:   Fitted PS model for unweighted combined sample                        ##
##   Formula.w: Fitted PS model for weighted combined sample. The default is NULL if  ##
##              the same as Formula (fitted to unweighted sample)                     ##
##   krnfun:    Kernel function                                                       ##
##   kw.mthd:   Choose PS model for PSAS method (can be vector).                      ##
##   rm.s:      Remove survey sample units not matched with any cohort units or not   ##
##              F: evenly distribute unmatched survey units weights to cohort units   ##
##              T: Remove unmatched survey units                                      ##
########################################################################################
simu_fun = function(samp.c, samp.s,
                    fm_fit.s, 
                    rm.s = F, 
                    fm_fit.o){
  rm.s <<- rm.s
  psa_dat = rbind(samp.c, samp.s)
  psa_dat$w   = c(rep(1, nrow(samp.c)), samp.s$wt)
  a = n_s/N
  psa_dat$w.s   = c(rep(1, nrow(samp.c)), samp.s$wt*a)
  
  wt_out = function(wt, a=1){
    ### Fit PS models and predict p(x)
    svyds = svydesign(ids =~1, strata=~trt, weight = wt, data = psa_dat)
    lgtreg = svyglm(fm_fit.s, family = binomial, design = svyds) 	
    beta = summary(lgtreg)$coeff[, 1]
    
    #Propensity scores 
    p.c = lgtreg$fitted.values[psa_dat$trt==1]
    p.s = lgtreg$fitted.values[psa_dat$trt==0]
    # Linear propensity scores
    p_score = lgtreg$linear.predictors
    # Linear propensity scores for the cohort
    p_score.c = p_score[psa_dat$trt==1]
    # Linear propensity scores for the survey sample
    p_score.s = p_score[psa_dat$trt==0]
    design.x = model.matrix(lgtreg)
    design.x.c = design.x[psa_dat$trt==1,]
    design.x.s = design.x[psa_dat$trt==0,]
    model.par = list(p.c=p.c, p.s=p.s, design.x.c=design.x.c, 
                     design.x.s=design.x.s)
    # IPSW
    ipsw = as.vector(1/exp(p_score.c))/a
    ipsw_beta = -ipsw*design.x.c
    #KW
    
    sgn_dist_mtx = outer(p_score.s, p_score.c, FUN = "-")
    # Bandwidths and weights and  # Kernal smoothing adjustment weights
    h = bw.nrd0(p_score.c)
    out = krnwt(sgn_dist_mtx, h, svy.wt = samp.s$wt, adj.mtx = 1, krnfun=dnorm, rm.s = F,
                w_beta=T, design.x.c=design.x.c, design.x.s=design.x.s)
    kw = out$psd.wt
    kw_beta = out$pw_beta
    colnames(kw_beta) = colnames(ipsw_beta)
    return(list(ipsw = ipsw, ipsw_beta = ipsw_beta,
                kw = kw, kw_beta = kw_beta,
                beta = beta, model.par = model.par))
  }
  # IPSW & KW.W
  wt_beta = wt_out(psa_dat$w)
  # IPSW
  ds.c = svydesign(ids=~1, data=samp.c, weights=wt_beta$ipsw, nest=TRUE)
  glm_out = svyglm(fm_fit.o, family=binomial, design = ds.c)
  beta_ipsw = summary(glm_out)$coeff[,1]
  mu_hat = glm_out$fitted.values
  design.x.o = model.matrix(glm_out)
  
  v_ipsw = v_Poisson(pw = wt_beta$ipsw, pw_beta = wt_beta$ipsw_beta, mu_hat = mu_hat, 
                     model.par = wt_beta$model.par,design.x.o = design.x.o, 
                     y = samp.c$y, svy.wt = samp.s$wt)
  # KW.W  
  ds.c = svydesign(ids=~1, data=samp.c, weights=wt_beta$kw, nest=TRUE)
  glm_out = svyglm(fm_fit.o, family=binomial, design = ds.c)
  beta_kw = summary(glm_out)$coeff[,1]
  mu_hat = glm_out$fitted.values
  v_kw = v_Poisson(pw = wt_beta$kw, pw_beta = wt_beta$kw_beta, mu_hat = mu_hat, 
                   model.par = wt_beta$model.par, design.x.o = design.x.o, 
                   y = samp.c$y, svy.wt = samp.s$wt)
  
  beta_est = rbind(beta_ipsw, beta_kw)
  beta_v = rbind(v_ipsw$v_gamma, v_kw$v_gamma)
  psd.wt = cbind(wt_beta$ipsw, wt_beta$kw)
  
  # IPSW.S & KW.S
  wt_beta = wt_out(psa_dat$w.s, a=a)
  # IPSW.S
  ds.c = svydesign(ids=~1, data=samp.c, weights=wt_beta$ipsw, nest=TRUE)
  glm_out = svyglm(fm_fit.o, family=binomial, design = ds.c)
  beta_ipsw = summary(glm_out)$coeff[,1]
  mu_hat = glm_out$fitted.values
  
  v_ipsw = v_Poisson(pw = wt_beta$ipsw, pw_beta = wt_beta$ipsw_beta, mu_hat = mu_hat, 
                     model.par = wt_beta$model.par, design.x.o = design.x.o, 
                     y = samp.c$y, svy.wt = samp.s$wt, a=a)
  # KW.S  
  ds.c = svydesign(ids=~1, data=samp.c, weights=wt_beta$kw, nest=TRUE)
  glm_out = svyglm(fm_fit.o, family=binomial, design = ds.c)
  beta_kw = summary(glm_out)$coeff[,1]
  mu_hat = glm_out$fitted.values
  v_kw = v_Poisson(pw = wt_beta$kw, pw_beta = wt_beta$kw_beta, mu_hat = mu_hat, 
                   model.par = wt_beta$model.par,design.x.o = design.x.o, 
                   y = samp.c$y, svy.wt = samp.s$wt, a=a)
  
  beta_est = rbind(beta_est, beta_ipsw, beta_kw)
  row.names(beta_est) = c("ipsw", "kw", "ipsw.s", "kw.s")
  beta_v = rbind(beta_v, v_ipsw$v_gamma, v_kw$v_gamma)
  row.names(beta_v) = c("ipsw", "kw", "ipsw.s", "kw.s")
  psd.wt = cbind(psd.wt, wt_beta$ipsw, wt_beta$kw)
  colnames(psd.wt) = c("ipsw", "kw", "ipsw.s", "kw.s")
  return(list(psd.wt = psd.wt, beta_est = beta_est, beta_v = beta_v))
}# end simu_fun



##############################################################################################
# FUNCION krnwt is an internal function of ksawt.                            # 
# INPUT                                                                                      #
#      sgn_dist_mtx
#      h
#      svy.wt
#      vary.h
kw.wt = function(p_score.c, p_score.s, h=NULL, svy.wt, krn, rm.s = F, Large = F,
                 w_beta=F, design.x.c=NULL, design.x.s=NULL){
  if(krn=="triang")h = bw.nrd0(p_score.c)/0.9*0.8586768
  if(krn=="dnorm"|krn=="dnorm_t")h = bw.nrd0(p_score.c)
  krnfun = get(krn)
  # create signed distance matrix    
  m = length(p_score.c)
  n = length(p_score.s)
  if(!Large){
    sgn_dist_mtx = outer(p_score.s, p_score.c, FUN = "-")
    kw_out = krnwt(sgn_dist_mtx = sgn_dist_mtx, h = h, svy.wt = svy.wt, 
                   krnfun = krnfun, rm.s = rm.s,
                   w_beta=w_beta, design.x.c=design.x.c, design.x.s=design.x.s)
    return(list(psd.wt  = kw_out$psd.wt, 
                pw_beta = kw_out$pw_beta, 
                sum_0.s = sum(kw_out$sum_0.s)))
  }else{
    psd.wt  = 0
    pw_beta = 0
    sum_0.s = 0
    grp_size =  floor(n/100)
    up = c(seq(0, n, grp_size)[2:100], n)
    lw = seq(1, n, grp_size)[-101]
    for(g in 1:100){
      sgn_dist_mtx = outer(p_score.s[lw[g]:up[g]], p_score.c, FUN = "-")
      if(!(is.null(design.x.s))){
        design.x.s.i = design.x.s[lw[g]:up[g],]
      }else{design.x.s.i=NULL}
      kw_out = krnwt(sgn_dist_mtx = sgn_dist_mtx, h = h, svy.wt = svy.wt[lw[g]:up[g]], 
                     krnfun = krnfun, rm.s = rm.s,
                     w_beta=w_beta, design.x.c=design.x.c, design.x.s=design.x.s.i)
      psd.wt  = psd.wt  + kw_out$psd.wt 
      pw_beta = pw_beta + kw_out$pw_beta 
      sum_0.s = sum_0.s + sum(kw_out$sum_0.s)      
    }
    return(list(psd.wt = psd.wt, 
                pw_beta = pw_beta, 
                sum_0.s = sum(sum_0.s)))
  }
}

krnwt = function(sgn_dist_mtx, h, svy.wt, krnfun, rm.s = F,
                 w_beta=F, design.x.c=NULL, design.x.s=NULL){
  m = ncol(sgn_dist_mtx)
  n = nrow(sgn_dist_mtx)
  krn_num = krnfun(sgn_dist_mtx/h)
  row.krn = rowSums(krn_num)
  sum_0.s = (row.krn==0)
  if(sum(sum_0.s)>0){
    warning('The input bandwidth h is too small. Please choose a larger one!')
    if(rm.s == T){
      warning(paste(sum(sum_0.s), "records in the prob sample were not used because of a small bandwidth"))
      row.krn[sum_0.s]=1
    }else{
      krn_num[sum_0.s,]= 1/m
      row.krn[sum_0.s] = 1
    }
  }
  
  krn = krn_num/row.krn
  # Final psuedo weights
  pswt_mtx = krn*svy.wt
  psd.wt = colSums(pswt_mtx)
  # Take derivative of beta
  pw_beta = NULL
  if(w_beta){
    if (is.null(design.x.c)|is.null(design.x.s)) stop("Need design matrix to claculate derivative of beta")
    n.beta = dim(design.x.c)[2]
    n_c =  dim(design.x.c)[1]
    pw_beta=matrix(0, n_c, n.beta)
    for(i in 2:n.beta){
      design.x.dist=outer(design.x.s[,i], design.x.c[,i], FUN="-")
      kij_beta = -krn_num*sgn_dist_mtx*design.x.dist/h/h
      row.kij_beta = rowSums(kij_beta)
      deriv1= (svy.wt/row.krn)%*%kij_beta
      deriv2= -(svy.wt*row.kij_beta/row.krn/row.krn)%*%krn_num
      pw_beta[,i] = deriv1+deriv2
    }
  }
  
  return(list(psd.wt = psd.wt, 
              pw_beta = pw_beta, 
              sum_0.s = sum(sum_0.s)))
  
} # end of krnwt



post_wt = function(sample, age_break, wt){
  sample$id = 1:nrow(sample)
  age_c = as.numeric(cut(sample$x1, breaks = age_break, include.lowest=T))
  sample$p_cell_s = age_c*100+sample$x2*10+ sample$x3_c
  sample = sample[order(sample $p_cell_s),]
  wt.sum = aggregate(wt, by = list(sample$p_cell_s), FUN=sum)
  age_c.p = as.numeric(cut(x1[y==1], breaks = age_break, include.lowest=T))
  p_cell_p = age_c.p*100+pop1$x2*10+pop1$x3_c
  pop.sum = aggregate(rep(1, N1), by = list(p_cell_p), FUN=sum)
  f = (pop.sum/wt.sum)[,2]
  sample$f = rep(f, as.numeric(table(sample$p_cell_s)))
  sample= sample[order(sample$id),]
  post.wt = wt*sample$f[order(sample$id)]
  post.wt
}

v_Poisson = function(pw, pw_beta, mu_hat, design.x.o, model.par, y, svy.wt, a=1){
  p.c = model.par$p.c
  p.s = model.par$p.s
  design.x.c = model.par$design.x.c
  design.x.s = model.par$design.x.s
  pi.c = p.c/(1-p.c)*a
  n_c = length(p.c); n_s = length(p.s)
  
  U_gamma = -t(pw*mu_hat*(1-mu_hat)*design.x.o)%*%design.x.o
  U_beta  = t((y-mu_hat)*design.x.o)%*%pw_beta
  S_gamma = matrix(0, nrow=ncol(design.x.c), ncol=ncol(design.x.o))
  S_beta  = -t(p.c*(1-p.c)*design.x.c)%*%design.x.c-  
    t(a*svy.wt* p.s*(1-p.s)*design.x.s)%*%design.x.s
  phi = rbind(cbind(U_gamma, U_beta),
              cbind(S_gamma, S_beta))
  phi_inv = solve(phi)
  
  U_gamma_inv = solve(U_gamma)
  S_beta_inv = solve(S_beta)
  b=-U_gamma_inv%*%U_beta%*%S_beta_inv
  phi_inv = rbind(cbind(U_gamma_inv, b),
                  cbind(S_gamma, S_beta_inv))
  # Calculate var(Phi)
  Phi_1  = as.matrix(cbind(pw*(y-mu_hat)*design.x.o, (1-p.c)*design.x.c)) 
  Phi_2  = as.matrix(cbind(matrix(0, nrow = nrow(design.x.s), ncol = ncol(design.x.o)), 
                           a*svy.wt* p.s*design.x.s))
  v1 = t(Phi_1)%*%(Phi_1*(1-c(pi.c)))
  v2 = t(Phi_2)%*%(Phi_2*(1-1/c(svy.wt)))
  v_Phi = v1+v2
  v_all = phi_inv%*%v_Phi%*%t(phi_inv)
  v_gamma = diag(v_all)[1:ncol(design.x.o)]
  return(list(v_mtx = v_all, v_gamma = v_gamma))  #, V_all = V_all  
} # End Sandwitch var

simu_fun_cox = function(pop, samp.c, samp.s,
                        fm_fit.s, 
                        rm.s = F, 
                        fm_fit.lg,
                        fm_fit.cox,
                        t_star){
  rm.s <<- rm.s
  psa_dat = rbind(samp.c, samp.s)
  psa_dat$w   = c(rep(1, nrow(samp.c)), samp.s$wt)
  a = n_s/N
  psa_dat$w.s   = c(rep(1, nrow(samp.c)), samp.s$wt*a)
  
  wt_out = function(wt, a=1){
    ### Fit PS models and predict p(x)
    svyds = svydesign(ids =~1, strata=~trt, weight = wt, data = psa_dat)
    lgtreg = svyglm(fm_fit.s, family = binomial, design = svyds) 	
    beta = summary(lgtreg)$coeff[, 1]
    
    #Propensity scores 
    p.c = lgtreg$fitted.values[psa_dat$trt==1]
    p.s = lgtreg$fitted.values[psa_dat$trt==0]
    # Linear propensity scores
    p_score = lgtreg$linear.predictors
    # Linear propensity scores for the cohort
    p_score.c = p_score[psa_dat$trt==1]
    # Linear propensity scores for the survey sample
    p_score.s = p_score[psa_dat$trt==0]
    design.x = model.matrix(lgtreg)
    design.x.c = design.x[psa_dat$trt==1,]
    design.x.s = design.x[psa_dat$trt==0,]
    model.par = list(p.c=p.c, p.s=p.s, design.x.c=design.x.c, 
                     design.x.s=design.x.s)
    # IPSW
    ipsw = as.vector(1/exp(p_score.c))/a
    ipsw_beta = -ipsw*design.x.c
    #KW
    sgn_dist_mtx = outer(p_score.s, p_score.c, FUN = "-")
    # Bandwidths and weights and  # Kernal smoothing adjustment weights
    h = bw.nrd0(p_score.c)
    out = krnwt(sgn_dist_mtx, h, svy.wt = samp.s$wt, adj.mtx = 1, krnfun=dnorm, rm.s = F,
                w_beta=F, design.x.c=design.x.c, design.x.s=design.x.s)
    kw = out$psd.wt
    #kw_beta = out$pw_beta
    #colnames(kw_beta) = colnames(ipsw_beta)
    return(list(ipsw = ipsw, ipsw_beta = ipsw_beta,
                kw = kw, #kw_beta = kw_beta,
                beta = beta, model.par = model.par))
  }
  est_prams = function(wt){
    ds.c = svydesign(ids=~1, data=samp.c, weights=as.formula(paste("~", wt)), nest=TRUE)
    prev =  svymean(~d, ds.c)[1]
    lg_coeff  = svyglm(fm_fit.lg, family=binomial, design = ds.c)$coeff
    cox_out = svycoxph(fm_fit.cox, design = ds.c)
    cox_coeff = cox_out$coefficients
    rel_hzd = exp(model.matrix(cox_out)%*%c(cox_coeff))
    
    cum_h = cum_bsln_hzd(dat = samp.c, wt = wt, 
                         rel_hzd = rel_hzd, 
                         t_star=t_star)
    gail_out = gail_h(pop = pop, samp = samp.c, wt = wt, 
                      rel_hzd = rel_hzd, 
                      t_star = t_star)
    
    return(list(prev       = prev, 
                lg_coeff   = lg_coeff,
                cox_coeff  = cox_coeff, 
                cum_h      = cum_h, 
                gail_1_ar  = gail_out$samp_1_ar,
                gail_cum_h = gail_out$cum_hzd))
  }
  # IPSW & KW.W
  wt_beta = wt_out(wt=psa_dat$w)
  # IPSW
  samp.c$ipsw = wt_beta$ipsw
  est_out    = est_prams(wt = "ipsw")
  prev_est   = est_out$prev
  lg_est     = est_out$lg_coeff
  cox_est    = est_out$cox_coeff
  cum_h_est  = est_out$cum_h
  gail_h_est = est_out$gail_cum_h
  gail_1_ar  = est_out$gail_1_ar
  # KW.W  
  samp.c$kww = wt_beta$kw
  est_out    = est_prams(wt = "kww")
  prev_est   = c(prev_est,       est_out$prev)
  lg_est     = rbind(lg_est,     est_out$lg_coeff)
  cox_est    = rbind(cox_est,    est_out$cox_coeff)
  cum_h_est  = rbind(cum_h_est,  est_out$cum_h)
  gail_h_est = rbind(gail_h_est, est_out$gail_cum_h)
  gail_1_ar  = rbind(gail_1_ar,  est_out$gail_1_ar)
  
  # IPSW.S & KW.S
  wt_beta = wt_out(psa_dat$w.s, a=a)
  # IPSW.S
  samp.c$ipsws = wt_beta$ipsw
  est_out    = est_prams("ipsws")
  prev_est   = c(prev_est,       est_out$prev)
  lg_est     = rbind(lg_est,     est_out$lg_coeff)
  cox_est    = rbind(cox_est,    est_out$cox_coeff)
  cum_h_est  = rbind(cum_h_est,  est_out$cum_h)
  gail_h_est = rbind(gail_h_est, est_out$gail_cum_h)
  gail_1_ar  = rbind(gail_1_ar,  est_out$gail_1_ar)
  # KW.S  
  samp.c$kws = wt_beta$kw
  est_out    = est_prams("kws")
  prev_est   = c(prev_est,       est_out$prev)
  lg_est     = rbind(lg_est,     est_out$lg_coeff)
  cox_est    = rbind(cox_est,    est_out$cox_coeff)
  cum_h_est  = rbind(cum_h_est,  est_out$cum_h)
  gail_h_est = rbind(gail_h_est, est_out$gail_cum_h)
  gail_1_ar  = rbind(gail_1_ar,  est_out$gail_1_ar)
  
  method = c("ipsw", "kw", "ipsw.s", "kw.s")
  row.names(lg_est)=row.names(cox_est)=row.names(cum_h_est)=method
  return(list(psd.wt     = samp.c[,(ncol(samp.c)-c(3:0))], 
              prev_est   = prev_est,
              lg_est     = lg_est, 
              cox_est    = cox_est,
              cum_h_est  = cum_h_est,
              gail_h_est = gail_h_est,
              gail_1_ar  = gail_1_ar))
}# end simu_fun_cox


samp.slct = function(seed, fnt.pop, n, Cluster=NULL, Clt.samp=NULL, dsgn, size = NULL, size.I = NULL){
  set.seed(seed)
  N = nrow(fnt.pop)
  size.Cluster = N/Cluster
  # one-ste sample design
  if(dsgn=="pps"){
    fnt.pop$x=size
    samp = sam.pps(fnt.pop,size, n)
  }
  # two-stage cluster sampling design (informative design at the second stage)
  if(dsgn == "srs-pps"){
    
    # calcualte the size for the second stage
    fnt.pop$x = size              
    
    #-- first stage: select clusters by srs
    # Clt.samp = 25
    # n = 1000
    index.psuI = sample(1:Cluster, Clt.samp, replac = F)
    index.psuI = sort(index.psuI)  #sort selected psus
    sample.I = fnt.pop[fnt.pop$psu %in% index.psuI,]
    sample.I$wt.I = Cluster/Clt.samp
    
    #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
    samp=NULL
    for (i in 1: Clt.samp){
      popn.psu.i= sample.I[sample.I$psu==index.psuI[i, 1],]
      size.II.i = sample.I[sample.I$psu==index.psuI[i, 1],"x"]
      samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  if(dsgn == "pps-pps"){
    
    # calcualte the size for the second stage
    fnt.pop$x = size               
    
    #-- first stage: select clusters by pps to size aggrated level p^a(or (1-p)^a)
    # Clt.samp = 25
    # n = 1000
    index.psuI = sam.pps(matrix(1:Cluster,,1),size.I, Clt.samp)
    index.psuI = index.psuI[order(index.psuI[,1]),]  #sort selected psus
    sample.I = fnt.pop[fnt.pop$psu %in% index.psuI[,1],]
    sample.I$wt.I = rep(index.psuI[,'wt'],each=size.Cluster)
    
    #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
    samp=NULL
    for (i in 1: Clt.samp){
      popn.psu.i= sample.I[sample.I$psu==index.psuI[i,1],]
      size.II.i = sample.I[sample.I$psu==index.psuI[i,1],"x"]
      samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  rownames(samp) = as.character(1:dim(samp)[1])
  return(samp)      
}


#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps<-function(popul,Msize, n){
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
    sam.pps.data=as.data.frame(popul[pps.samID,])
    names(sam.pps.data) = names(popul)
  }else{sam.pps.data=popul[pps.samID,]}
  sam.pps.data$wt=sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}

