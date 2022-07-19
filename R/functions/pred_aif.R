# get aif prediction based on gnm model results
# input model result and data after tmax
pred_aif <- function(mod_res,data){
    est_ac = coef(mod_res) # get coefficients from model results
    # calculate predicted aif using tri-exponential formula sum_i{exp(a_i+c_i*t})
    newtime = seq(0,max(data$time),0.01)
    pred_conc_plot = exp(est_ac[1]+est_ac[2]*newtime)+exp(est_ac[3]+est_ac[4]*newtime)+ exp(est_ac[5]+est_ac[6]*newtime)
    pred_conc = exp(est_ac[1]+est_ac[2]*data$time)+exp(est_ac[3]+est_ac[4]*data$time)+ exp(est_ac[5]+est_ac[6]*data$time)
    pred = list(plot = tibble(time_plot = newtime + unique(data$tmax), # for smooth plot
                              pred_plot = pred_conc_plot), 
                rsd = tibble(time_dsc = data$time, 
                             time = data$time +unique(data$tmax),
                             pred = pred_conc))
  return(pred)
}