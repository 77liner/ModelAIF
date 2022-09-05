# library(segmented)

#conc=data$count_asc[1] %>% as.data.frame()
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

# for dispersion aif
# find_t0_disp <- function(conc){
#   fit_lm = lm(disp_aif ~ time, data = conc)  # intercept-only model
#   fit_segmented = segmented(fit_lm, seg.Z = ~time, npsi = 1)  # Two change   points along x
#   t_0 = fit_segmented$psi[2]
#   conc = conc %>% mutate(t0 = t_0)
#   return(list(data = conc, segmod = fit_segmented, lmmod = fit_lm))
# }