
####
c=data$count[1] %>% as.data.frame()

c %>% mutate(tmax = time[which(aif == max(aif))]) ## tmax=0.666667



###########################
### ASE part
###########################
ase = c %>% 
  filter(time<=time[which(aif == max(aif))])

fit_lm = lm(aif ~ time, data = ase) # intercept-only model
plot(fit_lm)
plot(ase$time,ase$aif)
fit_segmented = segmented(fit_lm, seg.Z = ~time, npsi = 1)  # one change points along x
t_0 = fit_segmented$psi[2] # t0=0.4569598
# t0s are the same for aif and counts
# fit_lm_count = lm(count ~ time, data = conc)  # intercept-only model
# fit_segmented_count = segmented(fit_lm, seg.Z = ~time, npsi = 1)  # Two change   points along x
# t_0_count = fit_segmented$psi[2]
x = c(ase$time[which.max(ase$time[ase$time<t_0])], last(ase$time)) ## the max time below t_0 and the last time
y = c(ase$aif[which.min(abs(x[1]-ase$time))], last(ase$aif)) ## aif of the max time below t_0 and the last aif
xout = filter(ase,between(time, x[1],x[2]))$time
pred_asc = approx(x,y,xout = xout,method = "linear")
plot(pred_asc)

pred = tibble(time = pred_asc$x, pred = pred_asc$y)
ase = ase %>% mutate(# t0_count = t_0_count,
  t0 = t_0)

ase_data = ase
ase_segmod = fit_segmented
ase_lmmod = fit_lm
ase_pred = tibble(time = pred_asc$x, pred = pred_asc$y)



###############################
###  DSC part
###############################
dsc = c %>% 
  filter(time>=time[which(aif == max(aif))]) %>%
  mutate(t_G = t_G-time[which(aif == max(aif))],
         time = time-time[which(aif == max(aif))])

t = dsc$time # time since tmax
t_G = dsc$t_G # time point put in gamma count since injection time
y.sum = dsc$count
delta = dsc$delta# time in the gamma counter
vol = dsc$vol
pf = dsc$parentFraction
bpr = dsc$bpr
disp = rep(1,nrow(dsc))
props = c(1/6,1/2)
seed = 1

#################
### gnm_mod1
###################
###start_val = initial_para(t,y.sum, props)

#######
# initial_para
#######

startvalue <- vector()
x=t # time since tmax
y=y.sum # count

plot(x,y)
# first part
minx = x[1:(length(x)*props[1])] 
miny = y[1:(length(y)*props[1])] 
plot(minx, miny)
# second part
midx = x[(length(x)*props[1]):(length(x)*props[2])] 
midy = y[(length(y)*props[1]):(length(y)*props[2])] 
plot(midx,midy)
# third part
maxx = x[(length(x)*props[2]):length(x)] 
maxy = y[(length(y)*props[2]):length(y)] 
plot(maxx,maxy)
# plot(minx,miny,type = "o")
# plot(maxx,maxy,type = "o")

# fit linear model on the third part
strip1 = lm(log(maxy) ~ maxx)
plot(strip1)
coef(strip1)
startvalue[5:6] = c(coef(strip1)[1],coef(strip1)[2])

# calculating the residuals for the second part
# the model is y=exp(a+bx)
resid23 = midy - exp(startvalue[5]+startvalue[6]*midx)
# regress residuals on second part
strip2 = lm(log(abs(resid23)) ~ midx)
startvalue[3:4] = c(coef(strip2)[1],coef(strip2)[2])

# calculating the residuals for the first part
resid12 = miny - exp(startvalue[3]+startvalue[4]*minx)

# regress residuals on first part
strip3 = lm(log(abs(resid12)) ~ minx)
startvalue[1:2] = c(coef(strip3)[1],coef(strip3)[2])
startvalue

#################
### gnm_mod1
###################
start_val = startvalue
set.seed(seed)
calibration = 0.003
#Add decay correction, metabolism correction
fit_res = gnm(y.sum ~ -1 + offset(log(delta)) + offset(log(vol)) + offset(log(disp)) + offset(-log(2)/20.364*t_G)# decay correction
              + offset(-log(calibration)*rep(1,length(t))) + # calibration
                offset(-log(pf))+# metabolism uncorrection
                offset(log(bpr))# bpr
              + log_sum_exp(t),family = poisson(link = "log"),start = start_val)
quiet_gnm_mod1 <- quietly(gnm_mod1)

####################
### gnm_prop
####################
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





###################
count_asc = ase_data
asc_mod = ase_segmod
asc_pred = ase_pred

dsc_mod = try_res$result# save model fit data after tmax
##############


############
### pred_aif
############

# get prediction after tmax, contain interpolated aif
#dsc_pred = map2(dsc_mod,dsc,~pred_aif(.x,.y))

est_ac = coef(dsc_mod) # get coefficients from model results
# calculate predicted aif using tri-exponential formula sum_i{exp(a_i+c_i*t})
newtime = seq(0,max(dsc$time),0.01)
# model is y=sum(1-3) exp(ai+bi*t)
pred_conc_plot = exp(est_ac[1]+est_ac[2]*newtime)+exp(est_ac[3]+est_ac[4]*newtime)+ exp(est_ac[5]+est_ac[6]*newtime)
pred_conc = exp(est_ac[1]+est_ac[2]*dsc$time)+exp(est_ac[3]+est_ac[4]*dsc$time)+ exp(est_ac[5]+est_ac[6]*dsc$time)
plot(pred_conc_plot, newtime)
plot(pred_conc, dsc$time)
pred = list(plot = tibble(time_plot = newtime + unique(dsc$tmax), # for smooth plot
                          pred_plot = pred_conc_plot), 
            rsd = tibble(time = dsc$time +unique(dsc$tmax),
                         pred = pred_conc))

dsc_pred=pred

#############

pred = rbind(asc_pred,dsc_pred$rsd[-1,]) # combine predicted aif

###########
### add_residual
##########
# rsd = map2(c, pred, ~add_residual(.x,.y))
inter_time = intersect(c$time, pred$time)
count = c %>% 
  filter(time %in% inter_time)
pred1 = pred %>% 
  filter(time %in% inter_time)

rsd = tibble(time = inter_time, residual = (pred1$pred-count$aif))

##################
### Plot aif
##################
count=data[1,]
x1 =10 
x2= 90
legend = TRUE

# get subset data according to the range of time
data = count$count[[1]] %>% filter(time <= x2 & time >= x1)
asc_pred = count$asc_pred[[1]] %>% filter(time <= x2 & time >= x1) ## pred after t_0 
dsc_pred = count$dsc_pred[[1]]$plot %>% filter(time_plot <= x2 & time_plot >= x1) ## newtime pred
patient = count$patient

# plot aif and poisson fits
plot(data$time,data$aif, lty = 1,type="p",pch = 15,col = "grey",cex=0.5, main = paste0("patient:",patient),ylim = c(0,max(data$aif,asc_pred$pred,dsc_pred$pred_plot)))
lines(asc_pred$time, asc_pred$pred, col = rgb(red = 0, green = 1, blue = 0, alpha = 0.8),pch = 20,type = "l",lwd = 1)
lines(dsc_pred$time_plot, dsc_pred$pred_plot, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.8),pch = 20,type = "l",lwd = 1)
if(legend){
  legend(x = "topright", legend=c("AIF","Poisson"),col=c("grey","red"), lty=c(NA,1,1),pch=c(20,NA,NA), cex=0.8)}

###############
### plot residual
###############
data = data[1,]
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

