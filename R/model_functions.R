# fit non-linear poisson model and return the model results
kinfitr_gnm <- function(t, # sampling time 
                        t_G = NULL, # time that put in the gamma counter
                        y.sum,  # count or concentration 
                        delta, # measure time in gamma counter
                        vol, # volume
                        calibration, # calibration of units and volumes
                        decay_corrected = TRUE, # concentration decay corrected?
                        meta_corrected = TRUE, # concentration metabolism corrected?
                        pf = NULL, # fitted parent fractions. Later can be changed to a vector of parameters of the parent fraction function
                        bpr = NULL,
                        disp){
  # if input has been decay corrected
  if(decay_corrected == TRUE & !is.null(t_G)){
    # fit the model count ~ -1 + offset(log(delta)) + offset(log(vol)) + offset(-log(2)/20.364*t_G) + log_sum_exp(t) + offset(-log(0.003)*rep(1,length(t))) # calibration
    model_res = gnm_prop(t,t_G,y.sum, calibration=calibration, pf = rep(1,length(t)),bpr=rep(1,length(t)),disp,delta,vol)
  } else {
    # if decay corrected indicator is not
    if(decay_corrected == FALSE){
      stop("For concentration input, whether it was decay (metabolism) corredted or not should be indicated")
    }
    # if no input for t_G then deliver error
    if(is.null(t_G)){
      stop("Need time at which samples are put in the gamma")
    }
  }
  
  # and has been decay and metabolism corrected
  if(meta_corrected == TRUE & !is.null(pf)){
    # guo hill function
    # fit the model count ~ -1 + offset(log(delta)) + offset(log(vol)) 
    # + offset(-log(2)/20.364*t_G)# decay uncorrection
    # + offset(-log(0.003)*rep(1,length(t))) # calibration
    # + offset(-log(pf))# metabolism uncorrection
    # + log_sum_exp(t) # sum of exponential term
    model_res = gnm_prop(t,t_G,y.sum,calibration=calibration, pf=pf,bpr=rep(1,length(t)),disp,delta,vol)
    
    #### add other parent fraction function
  } else if (meta_corrected == TRUE){
    if(is.null(pf)){
      stop("Need fitted parent fraction")
    }
    if(meta_corrected == FALSE){
      stop("For concentration input, whether it was decay (metabolism) corredted or not should be indicated")
    }
  }
  return(model_res)
}

#conc=data$count_dsc[1] %>% as.data.frame()
# x = conc$time, y = conc$count, props = c(1/6,1/2)

initial_para <- function(x,y,props){
  startvalue <- vector()
  
  # first part
  minx = x[1:(length(x)*props[1])] 
  miny = y[1:(length(y)*props[1])] 
  # second part
  midx = x[(length(x)*props[1]):(length(x)*props[2])] 
  midy = y[(length(y)*props[1]):(length(y)*props[2])] 
  # third part
  maxx = x[(length(x)*props[2]):length(x)] 
  maxy = y[(length(y)*props[2]):length(y)] 
  
  # plot(minx,miny,type = "o")
  # plot(maxx,maxy,type = "o")
  
  # fit linear model on the third part
  strip1 = lm(log(maxy) ~ maxx)
  startvalue[5:6] = c(coef(strip1)[1],coef(strip1)[2])
  
  # calculating the residuals for the second part
  resid23 = midy - exp(startvalue[5]+startvalue[6]*midx)
  # regress residuals on second part
  strip2 = lm(log(abs(resid23)) ~ midx)
  startvalue[3:4] = c(coef(strip2)[1],coef(strip2)[2])
  
  # calculating the residuals for the first part
  resid12 = miny - exp(startvalue[3]+startvalue[4]*minx)
  
  # regress residuals on first part
  strip3 = lm(log(abs(resid12)) ~ minx)
  startvalue[1:2] = c(coef(strip3)[1],coef(strip3)[2])
  
  return(startvalue)
}

# nonlinear term definition for gnm() formula
log_sum_exp <- function(expression, inst = NULL){
  list(predictors = list(a1 = 1, c1 = 1, a2 = 1, c2 = 1, a3 = 1, c3 = 1),
       variables = list(substitute(expression)),
       term = function(predLabels, varLabels) {
         paste("log(exp(", predLabels[1],"+",predLabels[2],"*",varLabels[1], ")+
               exp(", predLabels[3],"+",predLabels[4],"*",varLabels[1], ")+
               exp(", predLabels[5],"+",predLabels[6],"*",varLabels[1], "))", sep = "")
       })
}
class(log_sum_exp) <- "nonlin"


# do gnm regression and return model results
gnm_mod1 <- function(t,t_G,y.sum,calibration,pf,bpr,disp,delta = delta,vol=vol, props = c(1/6,1/2),seed = 1){
  # find appropraite starting values
  start_val = initial_para(t,y.sum, props)
  set.seed(seed)
  #Add decay correction, metabolism correction
  fit_res = gnm(y.sum ~ -1 + offset(log(delta)) + 
                  offset(log(vol)) + 
                  offset(log(disp)) + 
                  offset(-log(2)/20.364*t_G) +# decay correction
                  offset(-log(calibration)) + # calibration
                  offset(-log(pf))+# metabolism uncorrection
                  offset(log(bpr)) +# blood to plasma ratio
                  log_sum_exp(t),family = poisson(link = "log"),start = start_val)
  return(fit_res)
  
}

# make codes keeping running when gnm regression failed to converge
quiet_gnm_mod1 <- quietly(gnm_mod1)


# sequentially find appropriate 
gnm_prop <- function(t,t_G,y.sum,calibration,pf = NULL,bpr = NULL,disp,delta,vol){
  try_res = quiet_gnm_mod1(t,t_G,y.sum,calibration=calibration, pf = pf,bpr=bpr,disp=disp,delta,vol)
  
  # if default (c(1/6,1/2)) not work, use cut c(1/60,1/10).
  props_new =c(1/60,1/10)
  while (length(try_res$messages)!=0 & props_new[1]<0.1){
    props_new[1] = props_new[1] + 0.02
    try_res = quiet_gnm_mod(t,y.sum,delta,vol,props = props_new)
  }
  
  if(length(try_res$messages)!=0){
    try_res = quiet_gnm_mod1(t,t_G,y.sum,calibration=calibration, pf=pf,bpr=bpr,disp=disp,delta,vol,props = c(1/90,6/90))
  }
  
  if(length(try_res$messages)!=0){
    try_res = quiet_gnm_mod1(t,t_G,y.sum,calibration=calibration, pf=pf,bpr=bpr,disp=disp,delta,vol,props = c(1/3,2/3))
  }
  
  
  # If still not working, try (m,1/2) where 1/6<m<1/3
  props_new =c(1/6,1/2)
  while (length(try_res$messages)!=0 & props_new[1]<1/3){
    props_new[1] = props_new[1]+0.01
    try_res = quiet_gnm_mod1(t,t_G,y.sum,calibration=calibration, pf=pf,bpr=bpr,disp=disp,delta,vol,props = props_new)
  }
  # If still not working, try (1/3,n) where 1/2<n<2/3
  props_new1 =c(1/3,1/2)
  while (length(try_res$messages)!=0 & props_new1[2]<2/3){
    props_new1[2] = props_new1[2]+0.01
    try_res = quiet_gnm_mod1(t,t_G,y.sum,calibration=calibration, pf=pf,bpr=bpr,disp=disp,delta,vol,props = props_new1)
  }
  return(try_res)
}

# find tmax
findtmax <- function (data){
  return(data %>% mutate(tmax = time[which(aif == max(aif))]))
}

# get data after tmax
slice_exp <- function(conc){
  conc = conc %>% 
    filter(time>=time[which(aif == max(aif))]) %>%
    mutate(t_G = t_G-time[which(aif == max(aif))],
           time = time-time[which(aif == max(aif))])
  return(conc)
}

# get data before tmax
slice_asc <- function(conc){
  conc = conc %>% 
    filter(time<=time[which(aif == max(aif))])
  return(conc)
}


# detect t0 and interpolate aif between t0 and tmax
# input data before tmax
# return data with t0 recorded, model result for finding t0, linear model result and predicted aif
acs_inter <- function(conc){
  fit_lm = lm(aif ~ time, data = conc)  # intercept-only model
  fit_segmented = segmented::segmented(fit_lm, seg.Z = ~time, npsi = 1)  # one change points along x
  t_0 = fit_segmented$psi[2] # t0
  # t0s are the same for aif and counts
  # fit_lm_count = lm(count ~ time, data = conc)  # intercept-only model
  # fit_segmented_count = segmented(fit_lm, seg.Z = ~time, npsi = 1)  # Two change   points along x
  # t_0_count = fit_segmented$psi[2]
  x = c(conc$time[which.max(conc$time[conc$time<t_0])], last(conc$time)) ## the max time below t_0 and the last time
  y = c(conc$aif[which.min(abs(x[1]-conc$time))], last(conc$aif)) ## aif of the max time below t_0 and the last aif
  xout = filter(conc,between(time, x[1],x[2]))$time
  pred_asc = approx(x,y,xout = xout,method = "linear")
  conc = conc %>% mutate(# t0_count = t_0_count,
    t0 = t_0)
  return(list(data = conc, 
              segmod = fit_segmented, 
              lmmod = fit_lm,
              pred = tibble(time = pred_asc$x, pred = pred_asc$y)))
}


# plot aif and poisson fits
scope_whole <- function(count,x1, x2,legend = TRUE){
  # get subset data according to the range of time
  data = count$count[[1]] %>% filter(time <= x2 & time >= x1)
  asc_pred = count$asc_pred[[1]] %>% filter(time <= x2 & time >= x1)
  dsc_pred = count$dsc_pred[[1]]$plot %>% filter(time_plot <= x2 & time_plot >= x1)
  patient = count$patient
  
  # plot aif and poisson fits
  plot(data$time,data$aif, lty = 1,type="p",pch = 15,col = "grey",cex=0.5, main = paste0("patient:",patient),ylim = c(0,max(data$aif,asc_pred$pred,dsc_pred$pred_plot)))
  lines(asc_pred$time, asc_pred$pred, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.8),pch = 20,type = "l",lwd = 1)
  lines(dsc_pred$time_plot, dsc_pred$pred_plot, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.8),pch = 20,type = "l",lwd = 1)
  if(legend){
    legend(x = "topright", legend=c("AIF","Poisson"),col=c("grey","red"), lty=c(NA,1,1),pch=c(20,NA,NA), cex=0.8)}
}

# plot residuals
plot_rsd <-function(data,x1, x2){
  patient = data$patient
  rsd = data$rsd[[1]]
  # get subset data according to the range of time
  rsd = rsd %>% 
    filter(time >= x1 & time <= x2) 
  # plot residual
  par(mfrow=c(1,2))
  plot(rsd$time,rsd$residual,pch = 20, cex = 1,col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),main = paste0("patient:",patient))
  abline(h=0,lty = 2,lwd = 0.8)
  
  # box plot of residual with median
  boxplot(rsd$residual)
  abline(h=median(rsd$residual),col = "red")
  text(1 - 0.4, median(rsd$residual), 
       labels = formatC(median(rsd$residual), format = "f", 
                        digits = 1),
       pos = 3, cex = 0.9, col = "red")
}
