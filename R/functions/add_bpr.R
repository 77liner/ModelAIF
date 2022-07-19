# calculate bpr
add_bpr <- function(count, conc){
  bprmodel = gam(bpr ~ s(time), data=count)
  bprpred = predict(bprmodel, newdata = list(time=conc$time))
  conc = conc %>% 
    mutate(bpr = bprpred,
           bpr = ifelse(Method == "Discrete", 1, bpr)) %>% 
    distinct(time,.keep_all = TRUE)
  return(conc)
}

# calculate aif
add_conc <- function(conc){
  conc = conc %>% 
    mutate(
      aif = count*exp(log(2)/20*time)*0.037/delta/vol*(0.07756609407608*exp(0.080728586340915*vol))/bpr*parentFraction,
      disp_aif = count*exp(log(2)/20*time)*0.037/delta/vol*(0.07756609407608*exp(0.080728586340915*vol))/bpr*parentFraction/disp_fct)
  return(conc)
}

# calculate residuals
add_residual <- function(count,pred){
  inter_time = intersect(count$time, pred$time)
  count = count %>% 
    filter(time %in% inter_time)
  pred = pred %>% 
    filter(time %in% inter_time)
  return(tibble(time = inter_time, residual = (pred$pred-count$aif)))
}


