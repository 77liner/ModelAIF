library(tidyverse)
library(scam)


# Readin data
filenames = list.files(path = "./data", pattern="*.csv") ##name of files
data = tibble(
  patient = substr(filenames,1,5), #the first 5 character of filenames
  count = map(filenames,~read.csv(paste0("./data/", .)) %>% select(-X))
)


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

# Seperate data
data = data %>% 
  group_by(patient) %>% 
  mutate(count = map(count,~findtmax(.x)), # the tmax is the one with max aif
         count_asc = map(count, ~slice_asc(.x)), # data before tmax
         count_dsc = map(count, ~slice_exp(.x))) # data after tmax and time=time-tmax; t_G=t_G-tmax

#c = data$count_dsc[1] %>% as.data.frame()

# non-parametric regression

non_regress = function(data = data, calibration = 0.003, disp = 1 ){
fit_res = scam(count ~ s(time,k =15, bs="mpd"),
               offset = log(delta)+log(vol)+log(disp)+(-log(2)/20.364*t_G)+(-log(calibration)*rep(1,length(time)))+(-log(parentFraction))+log(bpr)
              ,family = poisson(link = "log"),data = data)

return(fit_res)
}

# fit model

fit_data = data %>% 
  group_by(patient) %>% 
  mutate(
    dsc_mod = map(count_dsc,~non_regress(.x)),
    dsc_pred =map(dsc_mod,~exp(.x$linear.predictors)),
    dsc_res = map(dsc_mod,~.x$residuals)
    )

# plots

for (i in 1:nrow(fit_data)){
  # plots of prediction
  plot_data = fit_data$count_dsc[i] %>% as.data.frame()
  plot_pred = fit_data$dsc_pred[i] %>% as.data.frame()
  patient = fit_data$patient[i]
  names(plot_pred)[1] = "pred"
  plot(plot_data$time, plot_pred$pred, lty = 1,type="p",pch = 15,col = "grey",
       cex=0.5, main = paste0("patient:",patient))
  
  # plots of residuals
  res = fit_data$dsc_res[i]%>% as.data.frame()
  names(res)[1] = "residual"
  boxplot(res$residual)
  abline(h=median(res$residual),col = "red")
  text(1 - 0.4, median(res$residual), 
       labels = formatC(median(res$residual), format = "f", 
                        digits = 3),
       pos = 3, cex = 0.9, col = "red")
}


