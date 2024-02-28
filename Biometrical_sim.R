library(foreach)
library(doParallel)
library(parallel)
library(boot)
library(purrr)
library(survival)
library(survminer)
library(data.table)
library(survRM2)
library(pseudo)
library(mvtnorm)
library(coda)
library(numDeriv)
library(matrixStats)
library(betareg)
library(stats)
library( numDeriv)
library(simcausal)
registerDoParallel(6)

source("Paper1_func.R")


# ###############################
# #Please set your model parameters
# ###############################
# #pi model parameters
# beta_0_pi<- -1 #-1 for setting 1 and 0.5 for setting 2
# beta_1_pi<- c(1,2,-1.5) #c(1,2,-1.5) for setting 1 and c(0,0,0) for setting 2
# #mu model parameters
# alpha_0<- -2
# alpha_1<-c(1.2,2)

#one para in beta dist
nu <-3
#paraset
para_set<-c(alpha_0,alpha_1,beta_0_pi,beta_1_pi,nu)
#length of window
tau <-30
# #number of subjects
# nsubj<-1500

#initial parameter sets used in estimation
initial_para_est <- c(-1,1,1,-1,1,1,0,2)

##########################################
#calculate bias MSE... based on 1000 sims
##########################################
sim <- 1000
diff <-numeric()
diff_abs <- numeric()
MSE<-numeric()
Mean_tau_T_est_no_censoring <-numeric()
ASE <- numeric()

RMST_diff <-numeric()
RMST_diff_abs <-numeric()
RMST_MSE<-numeric()
RMST_predicted <- numeric()
RMST_ASE<-numeric()

diff_MI<-numeric()
diff_MI_abs <-numeric()
MSE_MI <- numeric()
MI_predicted<-numeric()
MI_ASE<-numeric()

diff_MI_indep<-numeric()
diff_MI_abs_indep <-numeric()
MSE_MI_indep <- numeric()
MI_predicted_indep<-numeric()
MI_ASE_indep<-numeric()

Pseudo_RMST_diff <- numeric()
Pseudo_RMST_diff_abs <- numeric()
Pseudo_RMST_MSE<-numeric()
Pseudo_predicted<-numeric()
Pseudo_ASE <- numeric()

Pseudo_RMST_diff_indep <- numeric()
Pseudo_RMST_diff_abs_indep <- numeric()
Pseudo_RMST_MSE_indep<-numeric()
Pseudo_predicted_indep<-numeric()
Pseudo_ASE_indep <- numeric()

empirical_var_restricted_mean<-numeric()
restricted_mean_CI_empirical <- numeric()
restricted_mean_CI_RMST_empirical<-numeric()
restricted_mean_CI_MI_empirical<-numeric()
restricted_mean_CI_MI_empirical_indep<-numeric()
restricted_mean_CI_RMST_pseudo_empirical <- numeric()
restricted_mean_CI_RMST_pseudo_empirical_indep <- numeric()

bias_coef<- list()
bias_sqr <- list()
empirical_se<-list()
cp <- list()

bias_coef_MI <- list()
se_MI <- list()
cp_MI<-list()
bias_coef_MI_indep <- list()
se_MI_indep <- list()
cp_MI_indep<-list()

converge_time<-numeric()
converge_time_indep<-numeric()
cens_num<-numeric()
cens_num_indep<-numeric()

bias_coef_EM<-list()
se_EM<-list()
cp_EM<-list()

bias_coef_EM_1<-list()
se_EM_1<-list()
cp_EM_1<-list()

diff_EM_1_indep<-numeric()
diff_EM_1_abs_indep<-numeric()
MSE_EM_1_indep<-numeric()
EM_1_ASE_indep<- numeric()
restricted_mean_CI_EM_1_empirical_indep<-numeric()

diff_EM_indep<-numeric()
diff_EM_abs_indep<-numeric()
MSE_EM_indep<-numeric()
EM_ASE_indep<- numeric()
restricted_mean_CI_EM_empirical_indep<-numeric()

Y_30<-numeric()

for (j in 1:sim){
  #simulate covariates
  print(j)
  set.seed(j)
  Z1<-runif(nsubj,0,1)
  Z2<-rbinom(n=nsubj, size=1, prob=0.7)
  Z3<-runif(nsubj,0,1)

  Z_mu<-cbind(Z1,Z2)
  Z_pi<-cbind(Z1,Z2,Z3)
  
  mu <- inv.logit(alpha_0+Z_mu%*%alpha_1)
  logit_mu <- logit(mu)
  alpha_b <- mu*nu
  beta_b <-(1-mu)*nu
  logit_pi <- beta_0_pi+Z_pi%*%beta_1_pi
  pi <- inv.logit(logit_pi)
  mu<-as.vector(mu)
  pi<-as.vector(pi)
  Mean_tau_T <- mu*tau*(1-pi)+tau*pi
  Y<-Generate_T()
  Y_original <-Y
  Y_30<-append(Y_30,sum(Y==30))

  ######################################
  #fit proposed model without censoring
  ######################################
  est<-nlminb(objective=nll_1,start = initial_para_est,control = list(iter.max=200),z_mu=Z_mu,z_pi=Z_pi,Y=Y)
  coef <- est$par
  coef_vcov<-Var_MI(coef,Y)
  se<-sqrt(diag(coef_vcov))
  empirical_se[[j]]<-se
  
  bias_coef[[j]] <- coef-para_set #bias of coefficients
  bias_sqr[[j]] <- bias_coef[[j]]^2 #squared bias
  
  #calculate cp of coefficients
  cp_num<-numeric()
  for (i in 1:length(para_set)){
    if (para_set[i]>(coef[i]-1.96*se[i]) & para_set[i]<(coef[i]+1.96*se[i])){cp_num[i]<-1}
    else{cp_num[i]<-0}
  }
  cp[[j]]<-cp_num
  
  mu_est <- inv.logit(coef[1]+Z_mu%*%coef[2:(ncol(Z_mu)+1)])
  pi_est <- inv.logit(coef[(ncol(Z_mu)+2)] + Z_pi%*%coef[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Mean_tau_T_est <- mu_est*tau*(1-pi_est)+tau*pi_est
  Mean_tau_T_est_no_censoring <- append(Mean_tau_T_est_no_censoring,Mean_tau_T_est)
  diff[j]<- sum(Mean_tau_T_est - Mean_tau_T)/length(Mean_tau_T) #bias
  diff_abs[j]<- sum(abs(Mean_tau_T_est - Mean_tau_T))/length(Mean_tau_T)
  MSE[j] <- sum((Mean_tau_T_est - Mean_tau_T)^2)/length(Mean_tau_T) #MSE
  
  
  #calculate variance of restricted mean and CI
  var_restricted_mean_est<-Var_restricted_mean_func(coef_vcov,coef)[[1]]
  ASE <- append(ASE,mean(sqrt(var_restricted_mean_est)))
  restricted_mean_lower <- Mean_tau_T_est-1.964739*sqrt(var_restricted_mean_est)
  restricted_mean_upper <- Mean_tau_T_est+1.964739*sqrt(var_restricted_mean_est)
  restricted_mean_CI <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>restricted_mean_lower[i] & Mean_tau_T[i]<restricted_mean_upper[i]){
      restricted_mean_CI <- restricted_mean_CI+1
    }
  }
  restricted_mean_CI_empirical[j] <- restricted_mean_CI/nsubj
  
  ##########################################
  #fit restricted mean model 
  ##########################################
  mean_tau_t <- numeric()
  for (m in 1:length(Y)){
    mean_tau_t[m] <- min(Y[m],tau)
  } #used in restricted mean model
  RMST_data <- data.frame(mean_tau_t,Z_pi) #Z_pi includes all covariates
  RMST <- lm(mean_tau_t ~ Z_pi, data=RMST_data) #RMST model
  RMST_coef <- RMST$coefficients
  #estimates of response (RMST model)
  RMST_mean_tau_t_est <- RMST_coef[1]+Z_pi%*%RMST_coef[2:(ncol(Z_pi)+1)]
  
  #estimates of the variance of restricted mean in RMST and CI of restricted mean
  RMST_prediction <- predict.lm(RMST,se.fit = TRUE,interval = "confidence")
  restricted_mean_CI_RMST <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>RMST_prediction$fit[i,2] & Mean_tau_T[i]<RMST_prediction$fit[i,3]){
      restricted_mean_CI_RMST <- restricted_mean_CI_RMST+1
    }
  }
  restricted_mean_CI_RMST_empirical[j] <- restricted_mean_CI_RMST/nsubj
  
  RMST_predicted <- append(RMST_predicted,RMST_prediction$fit[,1])
  RMST_diff[j] <- sum(RMST_mean_tau_t_est - Mean_tau_T)/length(Mean_tau_T) #bias
  RMST_diff_abs[j] <- sum(abs(RMST_mean_tau_t_est - Mean_tau_T))/length(Mean_tau_T) #bias
  RMST_MSE[j] <- sum((RMST_mean_tau_t_est - Mean_tau_T)^2)/length(Mean_tau_T)
  RMST_ASE <- append(RMST_ASE,mean(RMST_prediction$se.fit))
  
  
  ##############################################################
  # Generate dependently censored data and fit proposed model using MI and EM
  ##############################################################
  percent <- 0.36
  Data_cens <- generate_de_cens_dataset(Y)
  row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
  cens_num<- append(cens_num,length(row_MI))
  
  MI_result <- MI_impute()
  bias_coef_MI[[j]] <- MI_result[[1]]-para_set[-length(para_set)]
  se_MI[[j]] <- sqrt(diag(MI_result[[2]]))
  
  #calculate cp (MI)
  cp_num_MI<-numeric()
  for (i in 1:(length(para_set)-1)){
    if (para_set[i]>(MI_result[[1]][i]-1.96*se_MI[[j]][i]) & para_set[i]<(MI_result[[1]][i]+1.96*se_MI[[j]][i])){cp_num_MI[i]<-1}
    else{cp_num_MI[i]<-0}
  }
  cp_MI[[j]]<-cp_num_MI
  
  #calcute bias and EMSE of restricted mean
  mu_est_MI <- inv.logit(MI_result[[1]][1]+Z_mu%*%MI_result[[1]][2:(ncol(Z_mu)+1)])
  pi_est_MI <- inv.logit(MI_result[[1]][(ncol(Z_mu)+2)]+Z_pi%*%MI_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
  MI_predicted <- append(MI_predicted,Mean_tau_T_est_MI)
  diff_MI[j] <- sum(Mean_tau_T_est_MI - Mean_tau_T)/length(Mean_tau_T)
  diff_MI_abs[j] <- sum(abs(Mean_tau_T_est_MI - Mean_tau_T))/length(Mean_tau_T)
  MSE_MI[j] <- sum((Mean_tau_T_est_MI - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #calculate variance of restricted mean and CI (MI)
  var_restricted_mean_MI_est<-Var_restricted_mean_func(MI_result[[2]],MI_result[[1]])[[1]]
  MI_ASE<- append(MI_ASE,mean(sqrt(var_restricted_mean_MI_est)))
  restricted_mean_MI_lower <- Mean_tau_T_est_MI-1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_MI_upper <- Mean_tau_T_est_MI+1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_CI_MI <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>restricted_mean_MI_lower[i] & Mean_tau_T[i]<restricted_mean_MI_upper[i]){
      restricted_mean_CI_MI <- restricted_mean_CI_MI+1
    }
  }
  restricted_mean_CI_MI_empirical[j] <- restricted_mean_CI_MI/nsubj
  
  ################ EM results for dependent censoring
  EM_result_1<-MI_impute_converge()
  w_1 <- EM_result_1[[3]]
  coef_EM_1 <-EM_result_1[[1]]
  bias_coef_EM_1[[j]] <- coef_EM_1-para_set
  var_EM_1 <- EM_var(w_1,coef_EM_1)
  se_EM_1[[j]] <- sqrt(diag(var_EM_1))
  #se_EM_1[[j]] <-sqrt(diag(MI_impute_converge_bootstrap(Data_cens,Z_mu,Z_pi))) #bootstrap se
  
  #calculate cp (MI)
  cp_num_EM_1<-numeric()
  for (i in 1:(length(para_set)-1)){
    if (para_set[i]>(EM_result_1[[1]][i]-1.96*se_EM_1[[j]][i]) & para_set[i]<(EM_result_1[[1]][i]+1.96*se_EM_1[[j]][i])){cp_num_EM_1[i]<-1}
    else{cp_num_EM_1[i]<-0}
  }
  cp_EM_1[[j]]<-cp_num_EM_1
  
  #calcute bias and EMSE of restricted mean
  mu_est_MI <- inv.logit(EM_result_1[[1]][1]+Z_mu%*%EM_result_1[[1]][2:(ncol(Z_mu)+1)])
  pi_est_MI <- inv.logit(EM_result_1[[1]][(ncol(Z_mu)+2)]+Z_pi%*%EM_result_1[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
  diff_EM_1_indep[j] <- sum(Mean_tau_T_est_MI - Mean_tau_T)/length(Mean_tau_T)
  diff_EM_1_abs_indep[j] <- sum(abs(Mean_tau_T_est_MI - Mean_tau_T))/length(Mean_tau_T)
  MSE_EM_1_indep[j] <- sum((Mean_tau_T_est_MI - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #calculate variance of restricted mean and CI (MI)
  var_restricted_mean_MI_est<-Var_restricted_mean_func(var_EM_1,EM_result_1[[1]])[[1]]
  EM_1_ASE_indep<- append(EM_1_ASE_indep,mean(sqrt(var_restricted_mean_MI_est)))
  restricted_mean_MI_lower <- Mean_tau_T_est_MI-1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_MI_upper <- Mean_tau_T_est_MI+1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_CI_MI <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>restricted_mean_MI_lower[i] & Mean_tau_T[i]<restricted_mean_MI_upper[i]){
      restricted_mean_CI_MI <- restricted_mean_CI_MI+1
    }
  }
  restricted_mean_CI_EM_1_empirical_indep[j] <- restricted_mean_CI_MI/nsubj
  

  ##############################################################
  # Fitted restricted mean model for data with dependent censoring
  ##############################################################
  pseudo <- pseudomean(Data_cens$X,Data_cens$Delta)
  pseudo_fit <- lm(pseudo ~ Z_pi )
  pseudo_RMST_coef<-pseudo_fit$coefficients
  Pseudo_RMST_mean_tau_t_est <- pseudo_RMST_coef[1]+Z_pi%*%pseudo_RMST_coef[2:(ncol(Z_pi)+1)]
  Pseudo_predicted<-append(Pseudo_predicted,Pseudo_RMST_mean_tau_t_est)
  Pseudo_RMST_diff[j] <- sum(Pseudo_RMST_mean_tau_t_est - Mean_tau_T)/length(Mean_tau_T)
  Pseudo_RMST_diff_abs[j] <- sum(abs(Pseudo_RMST_mean_tau_t_est - Mean_tau_T))/length(Mean_tau_T)
  Pseudo_RMST_MSE[j] <- sum((Pseudo_RMST_mean_tau_t_est - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #estimates of the variance of restricted mean in RMST (when censoring is present)
  RMST_prediction_pseudo <- predict.lm(pseudo_fit,se.fit = TRUE,interval = "confidence")
  restricted_mean_CI_RMST_pseudo <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>RMST_prediction_pseudo$fit[i,2] & Mean_tau_T[i]<RMST_prediction_pseudo$fit[i,3]){
      restricted_mean_CI_RMST_pseudo <- restricted_mean_CI_RMST_pseudo+1
    }
  }
  Pseudo_ASE<- append(Pseudo_ASE,mean(RMST_prediction_pseudo$se.fit))
  restricted_mean_CI_RMST_pseudo_empirical[j] <- restricted_mean_CI_RMST_pseudo/nsubj
  
  
  ##############################################################################
  #Fit the proposed model for independently censored data using MI and EM
  ##############################################################################
  Data_cens <- generate_cens_dataset(Y)
  row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
  
  cens_num_indep<-numeric()
  cens_num_indep <- append(cens_num_indep,length(row_MI))
  MI_result <-MI_impute() #MI_impute()
  bias_coef_MI_indep[[j]] <- MI_result[[1]]-para_set[-length(para_set)]
  se_MI_indep[[j]] <- sqrt(diag(MI_result[[2]]))
  #calculate cp (MI)
  cp_num_MI_indep<-numeric()
  for (i in 1:(length(para_set)-1)){
    if (para_set[i]>(MI_result[[1]][i]-1.96*se_MI_indep[[j]][i]) & para_set[i]<(MI_result[[1]][i]+1.96*se_MI_indep[[j]][i])){cp_num_MI_indep[i]<-1}
    else{cp_num_MI_indep[i]<-0}
  }
  cp_MI_indep[[j]]<-cp_num_MI_indep
  
  #calcute bias and EMSE of restricted mean
  mu_est_MI <- inv.logit(MI_result[[1]][1]+Z_mu%*%MI_result[[1]][2:(ncol(Z_mu)+1)])
  pi_est_MI <- inv.logit(MI_result[[1]][(ncol(Z_mu)+2)]+Z_pi%*%MI_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
  MI_predicted <- append(MI_predicted,Mean_tau_T_est_MI)
  diff_MI_indep[j] <- sum(Mean_tau_T_est_MI - Mean_tau_T)/length(Mean_tau_T)
  diff_MI_abs_indep[j] <- sum(abs(Mean_tau_T_est_MI - Mean_tau_T))/length(Mean_tau_T)
  MSE_MI_indep[j] <- sum((Mean_tau_T_est_MI - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #calculate variance of restricted mean and CI (MI)
  var_restricted_mean_MI_est<-Var_restricted_mean_func(MI_result[[2]],MI_result[[1]])[[1]]
  MI_ASE_indep<- append(MI_ASE_indep,mean(sqrt(var_restricted_mean_MI_est)))
  restricted_mean_MI_lower <- Mean_tau_T_est_MI-1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_MI_upper <- Mean_tau_T_est_MI+1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_CI_MI <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>restricted_mean_MI_lower[i] & Mean_tau_T[i]<restricted_mean_MI_upper[i]){
      restricted_mean_CI_MI <- restricted_mean_CI_MI+1
    }
  }
  restricted_mean_CI_MI_empirical_indep[j] <- restricted_mean_CI_MI/nsubj
  
  
  
  ########################################
  #Fit our proposed model using EM algorithm
  EM_result<-MI_impute_converge()
  w <- EM_result[[3]]
  coef_EM <-EM_result[[1]]
  bias_coef_EM[[j]] <- coef_EM-para_set
  var_EM <- EM_var(w,coef_EM)
  se_EM[[j]] <- sqrt(diag(var_EM))
  cp_num_EM<-numeric()
  for (i in 1:(length(para_set)-1)){
    if (para_set[i]>(EM_result[[1]][i]-1.96*se_EM[[j]][i]) & para_set[i]<(EM_result[[1]][i]+1.96*se_EM[[j]][i])){cp_num_EM[i]<-1}
    else{cp_num_EM[i]<-0}
  }
  cp_EM[[j]]<-cp_num_EM
  
  #calcute bias and EMSE of restricted mean
  mu_est_MI <- inv.logit(EM_result[[1]][1]+Z_mu%*%EM_result[[1]][2:(ncol(Z_mu)+1)])
  pi_est_MI <- inv.logit(EM_result[[1]][(ncol(Z_mu)+2)]+Z_pi%*%EM_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
  diff_EM_indep[j] <- sum(Mean_tau_T_est_MI - Mean_tau_T)/length(Mean_tau_T)
  diff_EM_abs_indep[j] <- sum(abs(Mean_tau_T_est_MI - Mean_tau_T))/length(Mean_tau_T)
  MSE_EM_indep[j] <- sum((Mean_tau_T_est_MI - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #calculate variance of restricted mean and CI (MI)
  var_restricted_mean_MI_est<-Var_restricted_mean_func(var_EM,EM_result[[1]])[[1]]
  EM_ASE_indep<- append(EM_ASE_indep,mean(sqrt(var_restricted_mean_MI_est)))
  restricted_mean_MI_lower <- Mean_tau_T_est_MI-1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_MI_upper <- Mean_tau_T_est_MI+1.964739*sqrt(var_restricted_mean_MI_est)
  restricted_mean_CI_MI <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>restricted_mean_MI_lower[i] & Mean_tau_T[i]<restricted_mean_MI_upper[i]){
      restricted_mean_CI_MI <- restricted_mean_CI_MI+1
    }
  }
  restricted_mean_CI_EM_empirical_indep[j] <- restricted_mean_CI_MI/nsubj
  
  ##############################################################
  # Fitted restricted mean model for data with independent censoring
  ##############################################################
  pseudo <- pseudomean(Data_cens$X,Data_cens$Delta)
  pseudo_fit <- lm(pseudo ~ Z_pi )
  pseudo_RMST_coef<-pseudo_fit$coefficients
  Pseudo_RMST_mean_tau_t_est <- pseudo_RMST_coef[1]+Z_pi%*%pseudo_RMST_coef[2:(ncol(Z_pi)+1)]
  Pseudo_predicted<-append(Pseudo_predicted,Pseudo_RMST_mean_tau_t_est)
  Pseudo_RMST_diff_indep[j] <- sum(Pseudo_RMST_mean_tau_t_est - Mean_tau_T)/length(Mean_tau_T)
  Pseudo_RMST_diff_abs_indep[j] <- sum(abs(Pseudo_RMST_mean_tau_t_est - Mean_tau_T))/length(Mean_tau_T)
  Pseudo_RMST_MSE_indep[j] <- sum((Pseudo_RMST_mean_tau_t_est - Mean_tau_T)^2)/length(Mean_tau_T)
  
  #estimates of the variance of restricted mean in RMST (when censoring is present)
  RMST_prediction_pseudo <- predict.lm(pseudo_fit,se.fit = TRUE,interval = "confidence")
  restricted_mean_CI_RMST_pseudo <- 0
  for(i in 1:nrow(Z_mu)){
    if(Mean_tau_T[i]>RMST_prediction_pseudo$fit[i,2] & Mean_tau_T[i]<RMST_prediction_pseudo$fit[i,3]){
      restricted_mean_CI_RMST_pseudo <- restricted_mean_CI_RMST_pseudo+1
    }
  }
  Pseudo_ASE_indep<- append(Pseudo_ASE_indep,mean(RMST_prediction_pseudo$se.fit))
  restricted_mean_CI_RMST_pseudo_empirical_indep[j] <- restricted_mean_CI_RMST_pseudo/nsubj
}


bias_coef_no_censoring_result <- rowMeans(data.frame(bias_coef)) 
bias_coef_MI_indep_result <- rowMeans(data.frame(bias_coef_MI_indep)) 
bias_coef_MI_dep_result <- rowMeans(data.frame(bias_coef_MI)) 
bias_coef_EM_indep_result <- rowMeans(data.frame(bias_coef_EM)) 
bias_coef_EM_dep_result <- rowMeans(data.frame(bias_coef_EM_1)) 

ase_coef_no_censoring_result <- rowMeans(data.frame(empirical_se)) 
ase_coef_MI_indep_result <- rowMeans(data.frame(se_MI_indep)) 
ase_coef_MI_dep_result <- rowMeans(data.frame(se_MI)) 
ase_coef_EM_indep_result <- rowMeans(data.frame(se_EM)) 
ase_coef_EM_dep_result <- rowMeans(data.frame(se_EM_1)) 

esd_coef_no_censoring_result <- rowSds(as.matrix(data.frame(bias_coef)))
esd_coef_MI_indep_result <- rowSds(as.matrix(data.frame(bias_coef_MI_indep)))
esd_coef_MI_dep_result <- rowSds(as.matrix(data.frame(bias_coef_MI)))
esd_coef_EM_indep_result <- rowSds(as.matrix(data.frame(bias_coef_EM)))
esd_coef_EM_dep_result <- rowSds(as.matrix(data.frame(bias_coef_EM_1)))

cp_coef_no_censoring_result <- rowMeans(data.frame(cp)) 
cp_coef_MI_indep_result <- rowMeans(data.frame(cp_MI_indep)) 
cp_coef_MI_dep_result <- rowMeans(data.frame(cp_MI)) 
cp_coef_EM_indep_result <- rowMeans(data.frame(cp_EM)) 
cp_coef_EM_dep_result <- rowMeans(data.frame(cp_EM_1)) 

RMST_two_part_no_censoring_result <-c(mean(diff),mean(MSE),mean(ASE),mean(restricted_mean_CI_empirical))
RMST_traditional_no_censoring_result <-c(mean(RMST_diff),mean(RMST_MSE),mean(RMST_ASE),mean(restricted_mean_CI_RMST_empirical))
RMST_two_part_indep_censoring_result <-c(mean(diff_MI_indep),mean(MSE_MI_indep), mean(MI_ASE_indep), mean(restricted_mean_CI_MI_empirical_indep))
RMST_traditional_indep_censoring_result <-c(mean(Pseudo_RMST_diff_indep),mean(Pseudo_RMST_MSE_indep),mean(Pseudo_ASE_indep),mean(restricted_mean_CI_RMST_pseudo_empirical_indep))
RMST_EM_indep_censoring_result <-c(mean(diff_EM_indep),mean(MSE_EM_indep), mean(EM_ASE_indep), mean(restricted_mean_CI_EM_empirical_indep))
RMST_two_part_dep_censoring_result <-c(mean(diff_MI),mean(MSE_MI), mean(MI_ASE), mean(restricted_mean_CI_MI_empirical))
RMST_traditional_dep_censoring_result <-c(mean(Pseudo_RMST_diff),mean(Pseudo_RMST_MSE),mean(Pseudo_ASE),mean(restricted_mean_CI_RMST_pseudo_empirical))
RMST_EM_dep_censoring_result <-c(mean(diff_EM_1_indep),mean(MSE_EM_1_indep),mean(EM_1_ASE_indep),mean(restricted_mean_CI_EM_1_empirical_indep))

results_data_table1 <- cbind(c("alpha_0","alpha_1","alpha_2","beta_0","beta_1","beta_2","beta_3"),
                      bias_coef_no_censoring_result,
                      bias_coef_EM_indep_result,
                      bias_coef_MI_indep_result, 
                      bias_coef_EM_dep_result,
                      bias_coef_MI_dep_result,
                      ase_coef_no_censoring_result,
                      ase_coef_EM_indep_result,
                      ase_coef_MI_indep_result,
                      ase_coef_EM_dep_result,
                      ase_coef_MI_dep_result,
                      esd_coef_no_censoring_result,
                      esd_coef_EM_indep_result,
                      esd_coef_MI_indep_result,
                      esd_coef_EM_dep_result,
                      esd_coef_MI_dep_result,
                      cp_coef_no_censoring_result,
                      cp_coef_EM_indep_result,
                      cp_coef_MI_indep_result,
                      cp_coef_EM_dep_result,
                      cp_coef_MI_dep_result )

results_data_table2 <- cbind(RMST_two_part_no_censoring_result,
RMST_traditional_no_censoring_result,
RMST_two_part_indep_censoring_result,
RMST_EM_indep_censoring_result,
RMST_traditional_indep_censoring_result,
RMST_two_part_dep_censoring_result,
RMST_EM_dep_censoring_result,
RMST_traditional_dep_censoring_result)



# write.csv(results_data_table1[-nrow(results_data_table1),],file = "Table_1_results_500_beta0.csv")
# write.csv(results_data_table2,file = "Table_2_results_500_beta0.csv")


