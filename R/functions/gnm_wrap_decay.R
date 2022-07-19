#start_val = startvalue


source("functions/initial_para.R")
source("functions/log_sum_exp.R")
library(gnm)


# do gnm regression and return model results
gnm_mod1 <- function(t,t_G,y.sum,pf,bpr,disp,delta = delta,vol=vol, props = c(1/6,1/2),seed = 1,calibration){
  # find appropraite starting values
  start_val = initial_para(t,y.sum, props)
  set.seed(seed)
  #Add decay correction, metabolism correction
  fit_res = gnm(y.sum ~ -1 + offset(log(delta)) + offset(log(vol)) + offset(log(disp)) + offset(-log(2)/20.364*t_G)# decay correction
                + offset(-log(calibration)*rep(1,length(t))) + # calibration of the machine
                  offset(-log(pf))+# metabolism uncorrection
                  offset(log(bpr))# blood to plasma ratio
                  + log_sum_exp(t),family = poisson(link = "log"),start = start_val)
  return(fit_res)
  
}

# make codes keeping running when gnm regression failed to converge
quiet_gnm_mod1 <- quietly(gnm_mod1)


# sequentially find appropriate 
gnm_prop <- function(t,t_G,y.sum,pf = NULL,bpr = NULL,disp,delta,vol, calibration = 0.003){
  try_res = quiet_gnm_mod1(t,t_G,y.sum,pf = pf,bpr=bpr,disp=disp,delta,vol, calibration= calibration)
  
  # if default (c(1/6,1/2)) not work, use cut c(1/60,1/10).
  props_new =c(1/60,1/10)
  while (length(try_res$messages)!=0 & props_new[1]<0.1){
    props_new[1] = props_new[1] + 0.02
    try_res = quiet_gnm_mod(t,y.sum,delta,vol,props = props_new)
  }
  
  if(length(try_res$messages)!=0){
    try_res = quiet_gnm_mod1(t,t_G,y.sum,pf=pf,bpr=bpr,disp=disp,delta,vol,props = c(1/90,6/90), calibration = calibration)
  }
  
  if(length(try_res$messages)!=0){
    try_res = quiet_gnm_mod1(t,t_G,y.sum,pf=pf,bpr=bpr,disp=disp,delta,vol,props = c(1/3,2/3), calibration = calibration)
  }

  
  # If still not working, try (m,1/2) where 1/6<m<1/3
  props_new =c(1/6,1/2)
  while (length(try_res$messages)!=0 & props_new[1]<1/3){
    props_new[1] = props_new[1]+0.01
    try_res = quiet_gnm_mod1(t,t_G,y.sum,pf=pf,bpr=bpr,disp=disp,delta,vol,props = props_new, calibration= calibration)
  }
  # If still not working, try (1/3,n) where 1/2<n<2/3
  props_new1 =c(1/3,1/2)
  while (length(try_res$messages)!=0 & props_new1[2]<2/3){
    props_new1[2] = props_new1[2]+0.01
    try_res = quiet_gnm_mod1(t,t_G,y.sum,pf=pf,bpr=bpr,disp=disp,delta,vol,props = props_new1, calibration= calibration)
  }
  return(try_res)
}