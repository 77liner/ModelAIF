
# get aif prediction based on gnm model results
# input model result and data after tmax
pred_aif <- function(mod_res,data){
  est_ac = coef(mod_res) # get coefficients from model results
  # calculate predicted aif using tri-exponential formula sum_i{exp(a_i+c_i*t})
  newtime = seq(0,max(data$time),0.01)
  pred_conc_plot = exp(est_ac[1]+est_ac[2]*newtime) + 
    exp(est_ac[3]+est_ac[4]*newtime)+ exp(est_ac[5]+est_ac[6]*newtime)
  pred_conc = exp(est_ac[1]+est_ac[2]*data$time) + 
    exp(est_ac[3]+est_ac[4]*data$time) + exp(est_ac[5]+est_ac[6]*data$time)
  pred = list(plot = tibble(time_plot = newtime + unique(data$tmax), # for smooth plot
                            pred_plot = pred_conc_plot), 
              rsd = tibble(time_dsc = data$time, 
                           time = data$time +unique(data$tmax),
                           pred = pred_conc))
  return(pred)
}

# readdata <- function(data_type = NULL){
#   filenames = list.files(path = "./data", pattern="*.csv") ##name of files
#   data = tibble(
#     patient = substr(filenames,1,5), #the first 5 character of filenames
#     count = map(filenames,~read.csv(paste0("./data/", .)) %>% select(-X))
#   )
#   return(data)
# }


# plot aif and poisson fits
scope_whole <- function(count,x1, x2,legend = TRUE){
  # get subset data according to the range of time
  data = count$count[[1]] %>% filter(time <= x2 & time >= x1)
  asc_pred = count$asc_pred[[1]] %>% filter(time <= x2 & time >= x1)
  dsc_pred = count$dsc_pred[[1]]$plot %>% filter(time_plot <= x2 & time_plot >= x1)
  patient = count$patient
  
  # plot aif and poisson fits
  plot(data$time,data$aif, lty = 1,type="p",pch = 15,col = "grey",cex=0.5, 
       main = paste0("patient:",patient),
       ylim = c(0,max(data$aif,asc_pred$pred,dsc_pred$pred_plot)))
  lines(asc_pred$time, asc_pred$pred, col = rgb(red = 1, green = 0, 
                                                blue = 0, alpha = 0.8),
        pch = 20,type = "l",lwd = 1)
  lines(dsc_pred$time_plot, dsc_pred$pred_plot, col = rgb(red = 1, green = 0, 
                                                        blue = 0, alpha = 0.8),
        pch = 20,type = "l",lwd = 1)
  if(legend){
    legend(x = "topright", legend=c("AIF","Poisson"),col=c("grey","red"), 
           lty=c(NA,1,1),pch=c(20,NA,NA), cex=0.8)}
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
  plot(rsd$time,rsd$residual,pch = 20, cex = 1,col = rgb(red = 0, green = 0, 
                                                         blue = 1, alpha = 0.5),
       main = paste0("patient:",patient))
  abline(h=0,lty = 2,lwd = 0.8)
  
  # box plot of residual with median
  boxplot(rsd$residual)
  abline(h=median(rsd$residual),col = "red")
  text(1 - 0.4, median(rsd$residual), 
       labels = formatC(median(rsd$residual), format = "f", 
                        digits = 1),
       pos = 3, cex = 0.9, col = "red")
}
