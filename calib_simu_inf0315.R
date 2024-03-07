#install.packages(c("survival", "survminer"))
#rm(list=ls())
#setwd("/Volumes/wangl29/Calibration")
setwd("/home/wangl29/Calibration")
library(survival)
library(survminer)
library(survey)
library(matrixStats)
library(MASS)
library(data.table)
source("simu_fun.R")
source("taylor_deviate.R")
source("jk_fun.R")
#seed = read.table("seed.txt", header = T)
#seed1 = seed[,1]
#seed2 = seed[,2]

seed1 = seednumber1[,1]
seed2 = seednumber2[,1]
path="/home/wangl29/Calibration/inf_36k_mis/"
#seed=runif(100, min=0, max=10000)
#pop_beta = matrix(1, 100, 2)
#pop_lg   = matrix(1, 100, 3)

#for(simu in 1:100){
set.seed(8291.27)
N=300000
#N=100
n_c = 300; 
n_s = 600
fc = n_c/N;fc*100; fs = n_s/N;fs*100
m_jk=30; n_jk=60
x_mu = c(0, 0, 0)
sd_x = c(4, 2, 2)
#x_rho = matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), 3, 3)
x_rho = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
x_sigma = x_rho*outer(sd_x, sd_x)
x_mtrx = as.data.frame(mvrnorm(N, mu = x_mu, Sigma = x_sigma))
x1 = x_mtrx[,1]
x2 = x_mtrx[,2]
x3 = x_mtrx[,3]
#x1_c = as.numeric(cut(x1, breaks = quantile(x1, probs = seq(0,1,0.2)), include.lowest = T))
#x2_c = as.numeric(cut(x2, breaks = quantile(x2, probs = seq(0,1,0.2)), include.lowest = T))
x1_c = as.numeric(cut(x1, breaks = quantile(x1, probs = c(0,0.3,0.6,1)), include.lowest = T))
x2_c = as.numeric(cut(x2, breaks = quantile(x2, probs = c(0,0.3,0.6,1)), include.lowest = T))
table(x1_c)/N
table(x2_c)/N
beta = c(log(-log(.85)/15),0.2, 0.2, 0.3)
n_beta=length(beta)-1
gamma = exp(cbind(1, x1, x2, x3)%*%c(beta))
alpha = 1
t_i = rweibull(N, shape=1, scale=1/gamma)
t0 = runif(N, min=0, max=1)
c1 = 15-t0
c2 = rweibull(N, shape=1, scale=-15/log(0.9))
mean(c2<=c1)
t = apply(cbind(t_i, c1, c2), 1, FUN=min)
d = as.numeric(t_i<=apply(cbind(c1, c2), 1, FUN=min))
mean(d)*100

t_d_gen = function(beta_delta, sigma){
  t_delta = cbind(1, x1, x2, x1*x2)%*%c(beta_delta)+rnorm(N, 0.01)
  t_delta[t_delta<0]=0; range(t_delta)
  t_d = t_i+t_delta
  t_tilde = apply(cbind(t_d, c1, c2), 1, FUN=min)
  d_tilde = as.numeric(t_d<=apply(cbind(c1, c2), 1, FUN=min))
  #### Imputation
  pop = data.frame(id=1:N, x1, x1_c, x2, x2_c, x3, t, d, 
                   t_delta, t_tilde, d_tilde)
  #pop$imp_grp = pop$x2_c*10+pop$x1_c
  pop_m = pop[pop$d_tilde==1,]
  #t_delta_m = aggregate(pop_m$t_delta, by=list(pop_m$imp_grp), FUN=mean)
  #pop_m = pop_m[order(pop_m$imp_grp),]
  #pop_m$t_delta.imp = rep(t_delta_m[,2], table(pop_m$imp_grp))
  #pop_m = pop_m[order(pop_m$id),]
  ## Check results
  #tab_imp = table(pop_m$imp_grp, pop_m$t_delta.imp)
  #imp_tab = cbind(row.names(tab_imp),
  #                colnames(tab_imp)[match(rowSums(tab_imp),colSums(tab_imp))])
  #sum(as.numeric(imp_tab[,1])==t_delta_m[,1])
  #sum(round(as.numeric(imp_tab[,2]), 6)==round(t_delta_m[,2],6))
  pop_m$t_delta.imp=lm(t_delta~x1*x2, data=pop_m)$fitted.values
  pop_m$t_delta.imp[pop_m$t_delta.imp<0]=0
  pop$t.imp=pop$t
  pop$t.imp[pop$d_tilde==1] = pop_m$t_tilde - pop_m$t_delta.imp
  range(pop$t.imp[pop$d_tilde==1])
  pop$t.imp[pop$t.imp<0]=0
  return(list(t_d_pop = pop[,c(9:12)],
              t_delta.m = mean(t_delta),
              d_delta.m = mean(d_tilde)*100))
}
beta_delta1 = c(2,  0.01, 0.02, 0.01); sigma1 = 0.01
beta_delta2 = c(10, 0.2,  0.2,  0.1);  sigma2 = 0.2
beta_delta3 = c(10, 0,  0,  0);  sigma3 = 0.2
delta1 = t_d_gen(beta_delta=beta_delta1, sigma=sigma1)
delta2 = t_d_gen(beta_delta=beta_delta2, sigma=sigma2)
delta3 = t_d_gen(beta_delta=beta_delta3, sigma=sigma3)
rbind(c(delta1$t_delta.m,range(delta1$t_d_pop$t_delta)),
      c(delta2$t_delta.m,range(delta2$t_d_pop$t_delta)),
      c(delta3$t_delta.m,range(delta3$t_d_pop$t_delta)))
matrix(c(delta1$d_delta.m,
         delta2$d_delta.m,
         delta3$d_delta.m),,1)/100

pop = data.frame(id=1:N, x1, x1_c, x2, x2_c, x3, t, d, w=1,
                 delta1$t_d_pop, delta2$t_d_pop, delta3$t_d_pop)
names(pop)[-c(1:9)]=paste0(rep(names(delta1$t_d_pop),3), ".",rep(1:3,each=4))
pop0 = pop[pop$d==0,]; N0 = nrow(pop0)
apply(pop[,c("d", paste0("d_tilde.",1:3))],2,mean)

fm_fit.s   = as.formula("trt~x1+x2")
fm_fit.cox = as.formula("Surv(t, d) ~ x1+x2+x3")
cox_pop = coxph(fm_fit.cox, data = pop)
pop_beta = cox_pop$coefficients;pop_beta

ds.matx= model.matrix(lm(rep(1, N)~x1*d+x2*d))

cox_pop.imp1 = coxph(as.formula("Surv(t.imp.1, d_tilde.1) ~ x1+x2+x3"), data = pop)$coefficients
cox_pop.imp2 = coxph(as.formula("Surv(t.imp.2, d_tilde.2) ~ x1+x2+x3"), data = pop)$coefficients
cox_pop.imp3 = coxph(as.formula("Surv(t.imp.3, d_tilde.3) ~ x1+x2+x3"), data = pop)$coefficients
rbind(cox_pop.imp1, cox_pop.imp2, cox_pop.imp3)
t_star=c(1:15)
gail_pop_out = gail_h(pop, samp=pop, wt="w", rel_hzd=exp(cox_pop$linear.predictors), t_star=t_star)
cum_h_pop = gail_pop_out$cum_hzd
lambda_star = lambda_star.pop(pop)

# PPS
#fm_ps       = as.formula("trt~x1+x2*d")
#fm_ps.tilde = as.formula("trt~x1+x2*d_tilde")
#fm_ps = c("trt~x1+x2*d",paste0("trt~x1+x2*d_tilde.", 1:3))
fm_ps = "trt~x1+x2"
gamma_c = c(-0.15, -0.75, 0.1, 0, -0.2) #x1, d, x2, x1:d, x2:d
#fm_ps       = as.formula("trt~x1+x2")
#gamma_c = c(-0.15, 0, 0.1, 0, 0) #x1, d, x2, x1:d, x2:d
n_gamma = sum(gamma_c!=0)
#gamma_c.0_test = seq(-3.6,-3.5, 0.0001)
#gamma_c_test = cbind(gamma_c.0_test,
#                     matrix(rep(gamma_c, length(gamma_c.0_test)),
#                            length(gamma_c.0_test),,byrow=T)
#)
#pi.c_test = exp(ds.matx%*%t(gamma_c_test))
#gamma_c.0 = gamma_c.0_test[which.min(abs(rep(1,N)%*%pi.c_test-N*fc))]  
#gamma_c = c(gamma_c.0, gamma_c); gamma_c  			#search for the intercept with specified fc
odds.c = exp(ds.matx[,-1]%*%gamma_c);#sum(pi.c) #once find the intercept, create pi.c for popn
pi.c = odds.c*N*fc/sum(odds.c)
max(pi.c)/min(pi.c)
quantile(pi.c)

gamma_s = c(0.07, 0, 0.07, 0, 0)
#gamma_s.0_test = seq(-5,-4.92, 0.0001)
#gamma_s_test = cbind(gamma_s.0_test,
#                     matrix(rep(gamma_s, length(gamma_s.0_test)),
#                            length(gamma_s.0_test),,byrow=T)
#)

#pi.s_test = exp(ds.matx%*%t(gamma_s_test));sign(max(pi.s_test)*min(pi.s_test))
#gamma_s.0 = gamma_s.0_test[which.min(abs(apply(pi.s_test, 2, max)/apply(pi.s_test, 2, min)-20))]
#gamma_s = c(gamma_s.0, gamma_s); gamma_s
odds.s = exp(ds.matx[,-1]%*% gamma_s);
pi.s = odds.s*N*fc/sum(odds.s)
max(abs(pi.s))/min(abs(pi.s)) #poisson sampling with identified intercept
quantile(pi.s)

#x0 = as.matrix(as.numeric(quantile(x1, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))))
x0 = cbind(as.numeric(quantile(x1)), as.numeric(quantile(x2)), as.numeric(quantile(x3)))[c(2:4),]

#NSIMU=3
#mth.names = c("Naive", "Cht", "Svy", paste0("IPSW.", 0:3),  paste0(c("clib", "calib_ipsw", "calib_ipsw.t", "imp", "impcalib_ipsw", "impcalib_ipsw.t"), ".d", rep(c(1:3), each=6)) )
mth.names = c("Naive", "Cht", "Svy", "IPSW", paste0(c("clib", "calib_ipsw", "imp", "impcalib_ipsw"), ".d", rep(c(1:3), each=4)) )
n.mth=length(mth.names)
beta_est    = array(0, c(NSIMU, n_beta, n.mth))
beta_var    = array(0, c(NSIMU, n_beta, n.mth))
Lambda_est  = array(0, c(NSIMU, length(t_star), n.mth))
Lambda_var  = array(0, c(NSIMU, length(t_star), n.mth))
LambdaG_est = array(0, c(NSIMU, length(t_star), n.mth))
LambdaG_var = array(0, c(NSIMU, length(t_star), n.mth))
l.absR_est    = array(0, c(NSIMU, length(t_star), n.mth))
l.absR_var    = array(0, c(NSIMU, length(t_star), n.mth))
l.absRG_est   = array(0, c(NSIMU, length(t_star), n.mth))
l.absRG_var   = array(0, c(NSIMU, length(t_star), n.mth))
m.absR_est    = array(0, c(NSIMU, length(t_star), n.mth))
m.absR_var    = array(0, c(NSIMU, length(t_star), n.mth))
m.absRG_est   = array(0, c(NSIMU, length(t_star), n.mth))
m.absRG_var   = array(0, c(NSIMU, length(t_star), n.mth))
h.absR_est    = array(0, c(NSIMU, length(t_star), n.mth))
h.absR_var    = array(0, c(NSIMU, length(t_star), n.mth))
h.absRG_est   = array(0, c(NSIMU, length(t_star), n.mth))
h.absRG_var   = array(0, c(NSIMU, length(t_star), n.mth))
skip_simu = NULL
simu=1
pi.c=rep(fc, N)
start.time=date()
N<<-N; N1<<-sum(pop$d) 
fm_fit.cox<<-fm_fit.cox 
n_be<<-n_beta; mth.names<<-mth.names; t_star<<-t_star
fm_ps<<-fm_ps

for(simu in c(simu:NSIMU)){
  
  #set.seed(seed1[simu])
  #samp.c.indx = runif(N)<pi.c; n_c=sum(samp.c.indx)       # poisson sampling
  #samp.c = pop[samp.c.indx,]
  #samp.c$wt = 1/pi.c[samp.c.indx]
  #n_d[simu]=sum(samp.c$d)
  samp.c0 = samp.slct(seed     = seed1[simu], 
                      fnt.pop  = pop, 
                      n        = N*fc, 
                      #Cluster  = Cluster, 
                      #Clt.samp = psu_c, 
                      dsgn     = "pps", 
                      #size.I   = size.I_c, 
                      size     = odds.c)
  ################################################################################################
  samp.s0 = samp.slct(seed     = seed2[simu], 
                      fnt.pop  = pop, 
                      n        = N*fs, 
                      #Cluster  = Cluster, 
                      #Clt.samp = psu_c, 
                      dsgn     = "pps", 
                      #size.I   = size.I_c, 
                      size     = odds.s)
  rbind(apply(pop    [,c("d", "d_tilde.1","d_tilde.2", "d_tilde.3")],2,mean),
        apply(samp.c0[,c("d", "d_tilde.1","d_tilde.2", "d_tilde.3")],2,mean),
        apply(samp.s0[,c("d", "d_tilde.1","d_tilde.2", "d_tilde.3")],2,mean))
  
  ##################################################
  #start.t=date()
  jk_out = jk_fun(samp.c=samp.c0, samp.s=samp.s0, m_jk=m_jk, n_jk=n_jk, recal.wt=F, sampling="noninf")
  #end.t=date()
  #start.t;end.t
  jk_var = function(theta_jk){
    mean.est = sapply(1:n.mth, function(i) apply(theta_jk[,,i], 2,mean))
    sum_sq_theta = lapply(1:n.mth, function(i) (t(theta_jk[,,i]) - mean.est[,i])^2)
    var.est = sapply(1:n.mth, function(i) c(rep((m_jk-1)/m_jk, m_jk), rep((n_jk-1)/n_jk, n_jk))%*%t(sum_sq_theta[[i]]))
    return(list(mean.est=mean.est,var.est=var.est))
  }
  beta.inf = jk_var(jk_out$beta_est)
  beta_est[simu,,]     = beta.inf$mean.est;  beta_var[simu,,]    = beta.inf$var.est   
  Ldma.inf = jk_var(jk_out$Lambda_est)
  Lambda_est[simu,,]   = Ldma.inf$mean.est;  Lambda_var[simu,,]  = Ldma.inf$var.est 
  LdmaG.inf = jk_var(jk_out$LambdaG_est)
  LambdaG_est[simu,,]  = LdmaG.inf$mean.est; LambdaG_var[simu,,] = LdmaG.inf$var.est 
  absR.inf = jk_var(jk_out$l.absR_est);  absRG.inf = jk_var(jk_out$l.absRG_est)
  l.absR_est[simu,,] = absR.inf$mean.est;  l.absR_var[simu,,] = absR.inf$var.est 
  l.absRG_est[simu,,]= absRG.inf$mean.est; l.absRG_var[simu,,]= absRG.inf$var.est 
  absR.inf = jk_var(jk_out$m.absR_est);  absRG.inf = jk_var(jk_out$m.absRG_est)
  m.absR_est[simu,,] = absR.inf$mean.est;  m.absR_var[simu,,] = absR.inf$var.est 
  m.absRG_est[simu,,]= absRG.inf$mean.est; m.absRG_var[simu,,]= absRG.inf$var.est 
  absR.inf = jk_var(jk_out$h.absR_est);  absRG.inf = jk_var(jk_out$h.absRG_est)
  h.absR_est[simu,,] = absR.inf$mean.est;  h.absR_var[simu,,] = absR.inf$var.est 
  h.absRG_est[simu,,]= absRG.inf$mean.est; h.absRG_var[simu,,]= absRG.inf$var.est 
  
  print(simu)
}
#end.time=date()


pop_param = list(NSIMU_tot = NSIMU_tot,
                 NSIMU = NSIMU,
                 mth.names = mth.names,
                 x0 = x0,
                 n_x0 = nrow(x0),
                 beta = beta,
                 gamma_c = gamma_c,
                 pop_beta = pop_beta,
                 cum_h_pop = cum_h_pop,
                 n_t = length(cum_h_pop)#,
                 #Lambda_t = Lambda_t
                 
)


beta_est.mtx    = cbind(simu_id = k, matrix(beta_est   , NSIMU, n_beta *n.mth))
beta_var.mtx    = cbind(simu_id = k, matrix(beta_var   , NSIMU, n_beta *n.mth))
Lambda_est.mtx  = cbind(simu_id = k, matrix(Lambda_est , NSIMU, length(t_star)*n.mth))
Lambda_var.mtx  = cbind(simu_id = k, matrix(Lambda_var , NSIMU, length(t_star)*n.mth))
LambdaG_est.mtx = cbind(simu_id = k, matrix(LambdaG_est, NSIMU, length(t_star)*n.mth))
LambdaG_var.mtx = cbind(simu_id = k, matrix(LambdaG_var, NSIMU, length(t_star)*n.mth))
l.absR_est.mtx  = cbind(simu_id = k, matrix(l.absR_est , NSIMU, length(t_star)*n.mth))
l.absR_var.mtx  = cbind(simu_id = k, matrix(l.absR_var , NSIMU, length(t_star)*n.mth))
l.absRG_est.mtx = cbind(simu_id = k, matrix(l.absRG_est, NSIMU, length(t_star)*n.mth))
l.absRG_var.mtx = cbind(simu_id = k, matrix(l.absRG_var, NSIMU, length(t_star)*n.mth))
m.absR_est.mtx  = cbind(simu_id = k, matrix(m.absR_est , NSIMU, length(t_star)*n.mth))
m.absR_var.mtx  = cbind(simu_id = k, matrix(m.absR_var , NSIMU, length(t_star)*n.mth))
m.absRG_est.mtx = cbind(simu_id = k, matrix(m.absRG_est, NSIMU, length(t_star)*n.mth))
m.absRG_var.mtx = cbind(simu_id = k, matrix(m.absRG_var, NSIMU, length(t_star)*n.mth))
h.absR_est.mtx  = cbind(simu_id = k, matrix(h.absR_est , NSIMU, length(t_star)*n.mth))
h.absR_var.mtx  = cbind(simu_id = k, matrix(h.absR_var , NSIMU, length(t_star)*n.mth))
h.absRG_est.mtx = cbind(simu_id = k, matrix(h.absRG_est, NSIMU, length(t_star)*n.mth))
h.absRG_var.mtx = cbind(simu_id = k, matrix(h.absRG_var, NSIMU, length(t_star)*n.mth))

write.table(beta_est.mtx,    paste0(path, "beta_est."   ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(beta_var.mtx,    paste0(path, "beta_var."   ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(Lambda_est.mtx,  paste0(path, "Lambda_est." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(Lambda_var.mtx,  paste0(path, "Lambda_var." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(LambdaG_est.mtx, paste0(path, "LambdaG_est.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(LambdaG_var.mtx, paste0(path, "LambdaG_var.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(l.absR_est.mtx,  paste0(path, "l.absR_est." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(l.absR_var.mtx,  paste0(path, "l.absR_var." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(l.absRG_est.mtx, paste0(path, "l.absRG_est.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(l.absRG_var.mtx, paste0(path, "l.absRG_var.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(m.absR_est.mtx,  paste0(path, "m.absR_est." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(m.absR_var.mtx,  paste0(path, "m.absR_var." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(m.absRG_est.mtx, paste0(path, "m.absRG_est.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(m.absRG_var.mtx, paste0(path, "m.absRG_var.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(h.absR_est.mtx,  paste0(path, "h.absR_est." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(h.absR_var.mtx,  paste0(path, "h.absR_var." ,  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(h.absRG_est.mtx, paste0(path, "h.absRG_est.",  k, ".txt"), sep = ",", row.names = F, col.names = F)
write.table(h.absRG_var.mtx, paste0(path, "h.absRG_var.",  k, ".txt"), sep = ",", row.names = F, col.names = F)

#if(simu<NSIMU)simu=simu-1
#est.beta = t(sapply(1:n.mth, function(i) apply(beta_est[1:simu,,i],2,mean)))
#row.names(est.beta)=mth.names
##round(est.beta,3)
#relb.beta=t((t(est.beta)-beta[-1])/beta[-1]*100)
#var.beta = t(sapply(1:n.mth, function(i) apply(beta_est[1:simu,,i],2,var)))
#row.names(var.beta)=mth.names
#tl.beta = t(sapply(1:n.mth, function(i) apply(beta_var[1:simu,,i],2,mean)))
#row.names(tl.beta)=mth.names
#mse.beta = t((t(est.beta)-beta[-1])^2)+var.beta
#round(cbind(relb.beta, var.beta*1e3, mse.beta*1e3),2)
##########################################################################################
#est.Lambda = t(sapply(1:n.mth, function(i) apply(Lambda_est[1:simu,,i],2,mean)))
#row.names(est.Lambda)=mth.names
##round(est.beta,3)
#relb.Lambda=t((t(est.Lambda)-cum_h_pop)/cum_h_pop*100)
#var.Lambda = t(sapply(1:n.mth, function(i) apply(Lambda_est[1:simu,,i],2,var)))
#row.names(var.Lambda)=mth.names
#mse.Lambda = t((t(est.Lambda)-cum_h_pop)^2)+var.Lambda
#round(cbind(relb.Lambda[,c(1,8,15)], var.Lambda[,c(1,8,15)]*1e5, mse.Lambda[,c(1,8,15)]*1e5),2)
#
##G
#est.LambdaG = t(sapply(1:n.mth, function(i) apply(LambdaG_est[1:simu,,i],2,mean)))
#row.names(est.LambdaG)=mth.names
#relb.LambdaG=t((t(est.LambdaG)-cum_h_pop)/cum_h_pop*100)
#var.LambdaG = t(sapply(1:n.mth, function(i) apply(LambdaG_est[1:simu,,i],2,var)))
#row.names(var.LambdaG)=mth.names
#mse.LambdaG = t((t(est.LambdaG)-cum_h_pop)^2)+var.LambdaG
#round(cbind(relb.LambdaG[,c(1,8,15)], var.LambdaG[,c(1,8,15)]*1e5, mse.LambdaG[,c(1,8,15)]*1e5),2)
#
###################################################################################
#absR_pop = c(1-exp(-outer(exp(x0[,3]%*%c(pop_beta)),cum_h_pop, FUN="*")))
#est.absR = t(sapply(1:n.mth, function(i) apply(absR_est[1:simu,,i],2,mean)))
#row.names(est.absR)=mth.names
##round(est.beta,3)
#relb.absR=t((t(est.absR)-absR_pop)/absR_pop*100)
#var.absR = t(sapply(1:n.mth, function(i) apply(absR_est[1:simu,,i],2,var)))
#row.names(var.absR)=mth.names
#mse.absR = t((t(est.absR)-absR_pop)^2)+var.absR
#round(cbind(relb.absR[,c(1,8,15)], var.absR[,c(1,8,15)]*1e5, mse.absR[,c(1,8,15)]*1e5),2)
#
##G
#est.absRG = t(sapply(1:n.mth, function(i) apply(absRG_est[1:simu,,i],2,mean)))
#row.names(est.absRG)=mth.names
#relb.absRG=t((t(est.absRG)-absR_pop)/absR_pop*100)
#var.absRG = t(sapply(1:n.mth, function(i) apply(absRG_est[1:simu,,i],2,var)))
#row.names(var.absRG)=mth.names
#mse.absRG = t((t(est.absRG)-absR_pop)^2)+var.absRG
#round(cbind(relb.absRG[,c(1,8,15)], var.absRG[,c(1,8,15)]*1e5, mse.absRG[,c(1,8,15)]*1e5),2)
#
##
###