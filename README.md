
# Aims

Here I aim to demonstrate the use of our work to model the arterial input function for PET imaging using count data. Firstly, I will go through how we arrange the data ready for these models.  Next, I will demonstrate the application of each model to this data.


# Preparation

## Libraries

```{r}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(knitr)
library(cowplot)
library(mgcv)
library(ggforce)
library(gnm)
library(data.table)
library(scam)
library(DT)
library(ggplot2)
library(kableExtra)
```

We'll also load the model functions as well as some helper functions which we'll use along the way.

```{r}
source("model_functions.R")
source("helper_functions.R")
```


## Loading data

I have saved the prepared AIF data in the Rawdata folder. This dataset consists of 10 patients.

```{r}
# load data
filenames = list.files(path = "../Rawdata", pattern="*.csv") ##name of files
data = tibble(
  patient = substr(filenames,1,5), #the first 5 character of filenames
  count = map(filenames,~read.csv(paste0("../Rawdata/", .)) %>% select(-X))
) %>% 
  mutate(count = map(count, ~.x %>% 
                          mutate(Discrete=Method=="Discrete")))
```


# Data

## General Approach

Conventionally, modelling of the input function is performed by converting all the count data into radioactivity.  We also estimate values of the whole-blood-to-plasma ratio (BPR) using blood radioactivity and plasma radioactivity values, which are either modelled or interpolated from a few data points. We also estimate values of the plasma parent fraction using measured values, which are either interpolated or modelled.  

In this approach, we model the input function using the original count data, i.e. before converting the data to radioactivity. However, we still require estimates of the BPR and parent fraction over time.  For this reason, the data must firstly be processed in a conventional manner with radioactivity to derive estimates of the BPR and plasma parent fraction, as well as for assessing the effects of dispersion correction.  

In the data below, then, there will both be count data, as well as processed data generated using the counts data converted to radioactivity.  The strategy of this model is to make use of the counts data to allow a us to model the data with a more correct error distribution, but this requires first processing the data in a more conventional manner, and using the counts once again when applying the model itself.


## Example Data

We'll start out by looking at one of the individuals within the dataset.  The data is as follows:

```{r}
subject_1 = data$count[1]%>% 
  as.data.frame()

datatable(subject_1)
```

- **Time** : time of drawing blood sample (i.e. not of measurement of radioactivity).
- **Method** : there were two kinds of methods of drawing blood sample. One is automatically drawn and measured by an autosampler machine, named $Continuous$; the other is drawn manually, named $Discrete$.
- **Count** : the number of recorded counts from the gamma counter.
- **aif** : arterial input function estimated in the conventional manner.
- **backg** : backgroud. The measuring instrument needs background correction, and this is the estimated background radioactivity in count.
- **bpr**: blood-to-plasma ratio. We were interested in the counts in plasma, however, the machine could only measure counts in blood. By measuring both radioactivity in plasma and blood in manual samples, we model blood-to-plasma ratio along with time. Then the bpr was used to correct samples collected by machine.
- **ParentFraction** : The fraction of radioactivity originating from the parent compound, and not from radioactive daughter compounds (i.e. metabolites).
- **dis_frc**: used for dispersion correction. Dispersion correction is calculated as a fraction multiplier for each time point using the original and corrected AIF radioactivity values. In this way, we incorporate dispersion correction into the count model.
- **t_G** : time of measuring counts in plasma for manual samples/counts in blood for automatic samples (i.e. not of drawing the blood sample).
- **vol** : volume of the plasma sample.
- **delta** : duration for which the sample is in the gamma tube (i.e. longer duration results in more counts).
- **tmax** : the time when AIF reaches the peak.
- **disp_aif** : arterial input function after dispersion correction.


Here I will plot the AIF calculated in the conventional manner, after transforming arterial blood measurements to their respective arterial plasma estimates using the blood-to-plasma ratio, and transforming arterial plasma values to their respective metabolite-corrected estimates using the plasma parent fraction.


```{r}
ggplot(subject_1, aes( x = time, y = disp_aif))+
  geom_point(shape = "diamond", size = 2)+
  labs(x = "Time", y = "AIF",
       title = "AIF over time for patient jdcs1")+
  theme(axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 2, size = 15))+
  theme_light()
```



As we can see, the AIF rises sharply, and then slowly descends. We name the time when AIF reaches the peak as $tmax$. For AIF before the peak, we used simple linear regression to model it. We were mainly interested in AIF after the peak. We would use parametric and non-parametric ways to model data after $tmax$.  For this reason, we split the curve into ascent (asc) and descent (dsc).

```{r}
# find tmax and slice the curve
data = data %>% 
  group_by(patient) %>% 
  mutate(count = map(count,~findtmax(.x)), # the tmax is the one with max aif
         count_asc = map(count, ~slice_asc(.x)), # data before tmax
         count_dsc = map(count, ~slice_exp(.x))) # data after tmax and time=time-tmax; t_G=t_G-tmax
```



# Modelling


## Offsets

For regression of count data, we use offsets in the model to describe, for instance, if a sample spends longer in the counter, in which case more counts will be recorded.


(1) The way we calculate AIF: $disp\_aif = \frac{count\times exp^{\frac{log(2)\times time}{20}}\times 0.037\times (0.07756609407608\times exp^{0.080728586340915\times vol})\times parentFraction}{bpr\times disp\_fct\times delta\times vol}$ 

(2) The way we set offsets in the code: `offset =log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+(-log(0.003))+(-0.0807*vol)+(-log(parentFraction))+log(bpr)`


This is a bit complicated to follow. Let's break that down:

- $delta$ : time in the gamma counter. $Counts$ were divided by $delta$ when we calculated $AIF$.
- $vol$ : volume of the plasma sample. $Counts$ were divided by $vol$ when we calculated $AIF$.
- $disp\_fct$ : dispersion correction factor
- $-log(2)/20.364*t\_G$ : decay correction. The radio tracer in the blood or plasma would keep decaying. We need to correct it back to the time we want. The formula for decay correction : $A_t = A_0 * e ^{-\lambda t}$, where $A_0$ is the original activity count at time zero, $A_t$ is the activity at time $t$, $\lambda$ is the decay constant, and $t$ is the elapsed time.The decay constant is $\frac{ln(2)}{t_{1/2}}$ where $t_{1/2}$ is the half-life of the radioactive material of interest.
- $(-log(0.037*0.0776*\exp{(-0.0807*vol)})$ : calibration of the measured counts from the gamma counter : $Y = 0.037\times0.07756609407608\times exp^{0.080728586340915\times vol}$. The exponential volume calibration is used only for the manual discrete samples, and not the continuous samples. $0.037$ is used for calibrating units: 1 $picocurie$ $=$ 0.037 $Bq$.
- $parentFraction$ : the estimated metabolite free fraction at each time point.
- $bpr$ : the estimated blood-to-plasma ratio at each time point. This is set to 1 for all plasma measurements, and only used for the whole blood measurements which are corrected to their estimated plasma values.


We will add a column for the volumetric calibration which is equal to 1 for the continuous samples.

```{r}
data <- data %>% 
  mutate(count_dsc = map(count_dsc, ~.x %>% 
                           mutate(calibration = ifelse(Method=="Continuous",
                                                  yes = 0.037 * 0.0776 ,
                                                  no = 0.037 * 0.0776 * exp(0.0807*vol)))))
```

Then we will also calculate the AIF based on the counts and the offsets in order to be able to visualise the measured values.

```{r}
create_aif_values <- function(count_dsc) {
  
  suppressWarnings(
    count_dsc <- count_dsc %>% 
      mutate( total_offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+
                   (-log(calibration)) +
                   (-log(parentFraction))+log(bpr)) %>% 
      mutate(calc_aif = exp(log(count) - total_offset)) %>% 
      mutate(calc_aif = ifelse(is.nan(calc_aif), yes=0, no=calc_aif))
  )
  
}

data <- data %>% 
  mutate(count_dsc = map(count_dsc, create_aif_values))
```




## 1. Parametric (tri-exponential) Poisson Regression

```{r}
pare_data = data %>% 
  group_by(patient) %>% 
  mutate(asc_res = map(count_asc, ~acs_inter(.x)), # detect t0 and interpolate aif between t0 and tmax
         count_asc = map(asc_res, ~.x$data), # add t0 to the data
         asc_mod = map(asc_res, ~.x$segmod), # save model that detecting t0
         
         dsc_mod = map(count_dsc, # fit nonlinear poisson regression for descending part
                       ~kinfitr_gnm(t = .x$time, # time since tmax
                                    t_G = .x$t_G, # time point put in gamma count since injection time
                                    y.sum = .x$count, 
                                    delta = .x$delta, # time in the gamma counter
                                    vol = .x$vol,
                                    calibration = .x$calibration,
                                    pf = .x$parentFraction,
                                    bpr = .x$bpr,
                                    disp = rep(1,nrow(.x))
                                    )
                       ),
         dsc_mod = map(dsc_mod, ~.x$result), # save model fit data after tmax
         asc_pred = map(asc_res, ~.x$pred), # get prediction before tmax
         # get prediction after tmax, contain interpolated aif
         dsc_pred = map2(dsc_mod,count_dsc,~pred_aif(.x,.y))
         
         ) %>% 
  select(-asc_res)
```

- <span style="color:grey">Grey dots: the continuous data with blood sample collected by machine</span>
- <span style="color:red">Red dots: the discrete data with blood sample collected by experimenters</span>
- <span style="color:darkorange2">Orange line: the predicted value for parametric regression model</span>


```{r}
for (i in 1:nrow(pare_data)){
   patient = pare_data$patient[i] 

  data_line = pare_data$count_dsc[i] %>% as.data.frame()
  
  para_data = pare_data$dsc_pred[[i]]$rsd %>% as.data.frame()
  
  # plot of continuous data
  con_data = data_line %>% filter(Method == "Continuous")
  plot(con_data$time, log(con_data$calc_aif), col= "grey",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = paste0("Parametric Regression for patient:", patient))
  legend(50,5,c("Continuous data"), text.col = "grey",bty = "n")
  
  # plot of discrete data
  dis_data = data_line %>% filter(Method == "Discrete")
  par(new = TRUE)
  plot(dis_data$time, log(dis_data$calc_aif), col= "red",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = "")
  legend(50,4.5,c("Discrete data"), text.col = "red",bty = "n")
 
  # plots of parametric regression
  par(new = TRUE)
  plot(para_data$time_dsc, log(para_data$pred), type = "l", col = "darkorange2", lwd = 2,
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = ""
       )
  legend(50,4,c("Predicted value"), text.col = "darkorange2",bty = "n")
  
}
```



## 2. Non-parametric Poisson Regression

For the non-parametric poisson regression, we used $SCAM$ function to achieve monotone decreasing smooths.  We also added a covariate to account for the discrepancy observed between discrete and continuous data.  

```{r, message = FALSE, warning=FALSE}
poisson_regress = function(data = data, k_value = 15){
fit_res = scam(count ~ s(time,k = k_value, bs="mpd")+ Method, 
               offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+
                 (-log(calibration)) +
                 (-log(parentFraction))+log(bpr),
              family = poisson(link = "log"), data = data)
return(fit_res)
}
```

Following, we showed fitness for non-parametric poisson regression.

- <span style="color:grey">Grey dots: the continuous data with blood sample collected by machine</span>
- <span style="color:red">Red dots: the discrete data with blood sample collected by experimenters</span>
- <span style="color:darkgreen">Green line: the predicted value for non-parametric poisson regression model with index variable (Method). The dashed line shows the predicted values for the continuous data, while the solid line shows the predicted values for the discrete data. We consider the discrete data to represent the relevant data-generating process, and the continuous data is less reliable owing to some factor. In this study, it was speculated that the continuous sampler may have been placed too close to participants' bodies, resulting in extra counts´.</span>


```{r, message = FALSE, warning=FALSE}

for (i in 1:nrow(data)){
  patient = data$patient[i] 

  data_line = data$count_dsc[i] %>% as.data.frame()
  plot_line = data_line %>% poisson_regress()
  
  time = data_line$time
  fitted = log(plot_line$fitted.values)
  offset = plot_line$model$`(offset)`
  method = data_line$Method
  discrete = plot_line$coefficients[2]
  df = cbind(time, fitted, offset, method, discrete) %>% as.data.frame()
  
   test_df = df %>%  
    mutate(
      time = as.numeric(time),
     fitted = as.numeric(fitted),
      offset = as.numeric(offset),
     discrete = as.numeric(discrete),
    aif = fitted - offset
    )
   
   # create a data frame for prediction 
   pred_data = data_frame(
     time = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
     Method = "Discrete",
     delta = 10,
     vol = 1,
     disp_fct = 1,
     t_G = time+7,
     parentFraction = 0.8,
     bpr = 1
   )
   
   # predict 10 data points for discrete data at 0.1-1 time intervel 
   pred_df = predict(plot_line, pred_data) %>% as.data.frame()
   time_10 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
   pred_df = cbind(time_10, pred_df)
   
   names(pred_df)[1] = "time"
   names(pred_df)[2] = "aif"
   
   dis_df = test_df %>% filter(method == "Discrete") %>% select(time,aif)
   bind_dis_df = rbind(pred_df, dis_df) %>% mutate(time = as.numeric(time)) %>% arrange(time)
  
   con_df = test_df %>% filter(method == "Continuous")
  
  # plot of continuous data
  con_data = data_line %>% filter(Method == "Continuous")
  plot(con_data$time, log(con_data$calc_aif), col= "grey",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = paste0("Poisson Regression for patient:", patient))
  legend(50,5,c("Continuous data"), text.col = "grey",bty = "n")
  
  # plot of discrete data
  dis_data = data_line %>% filter(Method == "Discrete")
  par(new = TRUE)
  plot(dis_data$time, log(dis_data$calc_aif), col= "red",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = "")
  legend(50,4.5,c("Discrete data"), text.col = "red",bty = "n")
  # plots of non-parametric regression
  
    par(new = TRUE)
  plot(bind_dis_df$time,bind_dis_df$aif,
       type = "l", col = "darkgreen", lwd = 2,
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = ""
       )
  legend(50,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  
    par(new = TRUE)
  plot(con_df$time,con_df$aif,
       type = "l",lty = 2, col = "darkgreen", lwd = 2,
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = ""
       )
  legend(50,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  
}

```



### Problem of over-dispersion

Although the predicted value shows poisson regression model fits well, when we plotted the Q-Q plot, it revealed the problem of over-dispersion. 

- <span style="color:red">Red lines: the reference line</span>

- <span style="color:grey">Grey bands: the reference bands </span>

```{r, message = FALSE, warning=FALSE}
for (i in 1:nrow(data)){
   patient = data$patient[i] 

  data_line = data$count_dsc[i] %>% as.data.frame()
  plot_line = data_line %>% poisson_regress()
  
   qq.gam(plot_line,rep=500, s.rep = 10, level =1,pch=19,cex=.2, main = "", rl.col = 2)
 mtext(paste0("Q-Q plot for patient:", patient), side = 3, line = -1, outer = TRUE)
  
}
```


## 3. Non-parametric Model with Negative Binomial distribution

To resolve the problem of over-dispersion, we changed poisson distribution to negative binomial distribution. Because the $SCAM$ package does not support negative binomial distribution, we also changed the function to $GAM$ from the `mgcv` package.  However, this means that our non-parametric model is no longer shape-constrained.

```{r, message = FALSE, warning=FALSE}
gam_nb_regress = function(data = data,k_value = 15){
fit_res = gam(count ~ s(time,k = k_value)+Method, 
               offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+
                 (-log(calibration)) +
                 (-log(parentFraction))+log(bpr),
              family = nb(link = "log"), method="REML",data = data)
return(fit_res)
}
```

- <span style="color:grey">Grey dots: the continuous data with blood sample collected by machine</span>
- <span style="color:red">Red dots: the discrete data with blood sample collected by experimenters</span>
- <span style="color:darkgreen">Green line: the predicted value for non-parametric negative binomial regression model with index variable (Method)</span>


```{r, message = FALSE, warning=FALSE}

for (i in 1:nrow(data)){
   patient = data$patient[i] 

  data_line = data$count_dsc[i] %>% as.data.frame()
  plot_line = data_line %>% gam_nb_regress()

  
  time = data_line$time
  fitted = log(plot_line$fitted.values)
  offset = plot_line$offset
  method = data_line$Method

  df = cbind(time, fitted, offset, method) %>% as.data.frame()
  
   test_df = df %>%  
    mutate(
      time = as.numeric(time),
     fitted = as.numeric(fitted),
      offset = as.numeric(offset),
    aif = fitted - offset
    )
   
   # create a data frame for prediction
   pred_data = data_frame(
     time = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
     Method = "Discrete",
     delta = 10,
     vol = 1,
     disp_fct = 1,
     t_G = time+7,
     parentFraction = 0.8,
     bpr = 1
   )
   
   # predict 10 data points for discrete data at 0.1-1 time intervel 
   pred_df = predict(plot_line, pred_data) %>% as.data.frame()
   time_10 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
   pred_df = cbind(time_10, pred_df)
   
   names(pred_df)[1] = "time"
   names(pred_df)[2] = "aif"
   
   dis_df = test_df %>% filter(method == "Discrete") %>% select(time,aif)
   bind_dis_df = rbind(pred_df, dis_df) %>% mutate(time = as.numeric(time)) %>% arrange(time)
  
   con_df = test_df %>% filter(method == "Continuous")
  
  # plot of continuous data
  con_data = data_line %>% filter(Method == "Continuous")
  plot(con_data$time, log(con_data$disp_aif), col= "grey",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = paste0("Regression for patient:", patient))
  legend(50,5,c("Continuous data"), text.col = "grey",bty = "n")
  
  # plot of discrete data
  dis_data = data_line %>% filter(Method == "Discrete")
  par(new = TRUE)
  plot(dis_data$time, log(dis_data$disp_aif), col= "red",
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = "")
  legend(50,4.5,c("Discrete data"), text.col = "red",bty = "n")
  # plots of non-parametric regression
  
    par(new = TRUE)
  plot(bind_dis_df$time,bind_dis_df$aif,
       type = "l", col = "darkgreen", lwd = 2,
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = ""
       )
  legend(50,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  
    par(new = TRUE)
  plot(con_df$time,con_df$aif,
       type = "l",lty = 2, col = "darkgreen", lwd = 2,
       xlab = "time(min)",
       ylab = "log(AIF)",
       xlim = c(0,90),
       ylim = c(-2,5),
       main = ""
       )
  legend(50,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  

}

```
### Solved problem of over-dispersion

- <span style="color:red">Red lines: the reference line</span>

- <span style="color:grey">Grey bands: the reference bands </span>

```{r, message = FALSE, warning=FALSE}
for (i in 1:nrow(data)){
   patient = data$patient[i] 

  data_line = data$count_dsc[i] %>% as.data.frame()
  plot_line = data_line %>% gam_nb_regress()
  
  qq.gam(plot_line,rep=500)
 mtext(paste0("Q-Q plot for patient:", patient), side = 3, line = -1, outer = TRUE)
}
```

### Problem of excessive wiggliness

While this model appears to have resolved the issue of the over-dispersion, we were dissatisfied with the fits, as they seem to be overfitting, with upward deviations along what we would expect to a practically monotonic descent.

### Solved problem of excessive wiggliness

To this end, we log-transformed time within the model in order to distribute the basis functions more evenly across the curve, with respect to both the number of measurements as well as the expected wiggliness (i.e. second derivative) of the curve.  First, we added back tmax for all patients' data. Then, we log-transformed time value and fit the model. Following shows the regression model, residual plot, and how the model fitted. 


<!-- ## A model across all patients -->

<!-- ```{r} -->
<!-- #Create a dataset for all patients' data -->

<!-- total_data = data_frame() -->

<!-- for(i in 1:nrow(data)){ -->
<!-- patient = data$patient[i]  -->
<!-- data_line = data$count_dsc[i] %>% as.data.frame() -->
<!-- p_data = cbind(patient, data_line) -->
<!-- total_data = rbind(total_data, p_data) -->
<!-- } -->

<!-- # dataset for manual data only -->
<!-- discrete_data = total_data %>%mutate(patient = as.factor(patient))%>% filter(Method == "Discrete") -->

<!-- ``` -->

<!-- ```{r} -->
<!-- gam_total_nb_regress = function(data = data){ -->
<!-- fit_res = gam(count ~ Method*patient+patient+s(time,k = 20)+s(time, patient, k=5, bs="fs"), -->
<!--                offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+(-log(0.003))+(-0.0807*vol) -->
<!--                +(-log(parentFraction))+log(bpr) -->
<!--               ,family = nb(link = "log"),data = data,method="REML") -->
<!-- return(fit_res) -->
<!-- } -->

<!-- total_data = total_data %>% mutate( -->
<!--   patient = as.factor(patient) -->
<!--   ) -->
<!-- ``` -->


<!-- ```{r} -->
<!-- total_line = total_data%>% gam_total_nb_regress() -->

<!-- plot(total_line, main = "Regression for all patients") -->
<!-- ``` -->

<!-- ### Q-Q plot for the model across all patients -->

<!-- ```{r} -->
<!--  qq.gam(total_line,rep=500) -->
<!--  mtext(paste0("Q-Q plot for patient:", patient), side = 3, line = -1, outer = TRUE) -->
<!-- ``` -->

<!-- ## Increase Smoothness -->

<!-- By using negative binomial regression, we dealed with the problem of over-dispersion and got model with good fit. However those models were too wiggly. We wanted to increase the smoothness by log-transform $time$. -->





- <span style="color:grey">Grey dots: the continuous data with blood sample collected by machine</span>
- <span style="color:red">Red dots: the discrete data with blood sample collected by experimenters</span>
- <span style="color:darkgreen">Green line: the predicted value for non-parametric negative binomial regression model with index variable (Method)</span>

```{r, message = FALSE, warning=FALSE}
tmax = c() 
nrow = c()
patients = c()
for(i in 1:nrow(data)){
   patient = data$patient[i] 

  patient_data = data$count_dsc[i]%>% as.data.frame()%>% mutate(
                      tmax = as.numeric(tmax),
                      time = as.numeric(time)
                      )
  patient_tmax = mean(patient_data$tmax)
  number = nrow(patient_data)
  patients = rbind(patients,patient)
  tmax = rbind(tmax, patient_tmax)
  nrow = rbind(nrow,number)
  # transform time : time =ln(time + tmax)
  time_data = patient_data %>% mutate(time = log(time+tmax))
  
  plot_line = time_data %>% gam_nb_regress()
  
  # plot of the model
  par(mfrow = c(1,1))
  
  # plot of predicted value
  time = time_data$time
  fitted = log(plot_line$fitted.values)
  offset = plot_line$offset
  method = time_data$Method

  df = cbind(time, fitted, offset, method) %>% as.data.frame() %>% mutate(
    fitted = as.numeric(fitted),
    offset = as.numeric(offset),
    aif = fitted-offset
  )
  con_df = df %>% filter(method == "Continuous")
  dis_df = df %>% filter(method == "Discrete")
  
  par(mfrow = c(1,1))
  con_data = time_data %>% filter(Method == "Continuous")
  plot(con_data$time, log(con_data$calc_aif), col= "grey",
       xlab = "log(time)",
       ylab = "log(AIF)",
       xlim = c(-1,5),
       ylim = c(-2,5),
       main = paste0("Regression for patient:", patient))
  legend(3,5,c("Continuous data"), text.col = "grey",bty = "n")
  
  # plot of discrete data
  dis_data = time_data %>% filter(Method == "Discrete")
  par(new = TRUE)
  plot(dis_data$time, log(dis_data$calc_aif), col= "red",
       xlab = "log(time)",
       ylab = "log(AIF)",
       xlim = c(-1,5),
       ylim = c(-2,5),
       main = "")
  legend(3,4.5,c("Discrete data"), text.col = "red",bty = "n")
  # plots of non-parametric regression
  
    par(new = TRUE)
  plot(con_df$time,con_df$aif,
       type = "l",lty = 2, col = "darkgreen", lwd = 2,
       xlab = "log(time)",
       ylab = "log(AIF)",
       xlim = c(-1,5),
       ylim = c(-2,5),
       main = ""
       )
  legend(3,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  
    par(new = TRUE)
  plot(dis_df$time,dis_df$aif,
       type = "l", col = "darkgreen", lwd = 2,
       xlab = "log(time)",
       ylab = "log(AIF)",
       xlim = c(-1,5),
       ylim = c(-2,5),
       main = ""
       )
  legend(3,4,c("Predicted value"), text.col = "darkgreen",bty = "n")
  

}

table = cbind(patients, tmax, nrow) %>% as.tibble()
colnames(table) = c("Patient", "Tmax", "Number_of_data")
table%>%knitr::kable(
  caption = "Tmax for each patient",
  align = "lcc"
)%>%
  kable_classic_2()
```

These curves look much more like what we would expect them to look like!  


And the QQ plots still appear to look good, shown below.

```{r, message = FALSE, warning=FALSE}
tmax = c() 
nrow = c()
patients = c()
for(i in 1:nrow(data)){
   patient = data$patient[i] 

  patient_data = data$count_dsc[i]%>% as.data.frame()%>% mutate(
                      tmax = as.numeric(tmax),
                      time = as.numeric(time)
                      )
  patient_tmax = mean(patient_data$tmax)
  number = nrow(patient_data)
  patients = rbind(patients,patient)
  tmax = rbind(tmax, patient_tmax)
  nrow = rbind(nrow,number)
  # transform time : time =ln(time + tmax)
  time_data = patient_data %>% mutate(time = log(time+tmax))
  
  plot_line = time_data %>% gam_nb_regress()
  
  # plot of the model
  par(mfrow = c(1,1))
  
   # Q-Q plot
  qq.gam(plot_line,rep=500)
  mtext(paste0("Q-Q plot for patient:", patient), side = 3, line = -1, outer = TRUE)

  

}
```


### Figure

```{r, fig.height=5, fig.width=6, message = FALSE, warning=FALSE}
create_preds <- function(data, patient_id) {
  
  patdat <- data %>% 
    filter(patient==patient_id) %>% 
    select(patient, count_dsc) %>% 
    unnest(count_dsc) %>% 
    mutate(tmax = as.numeric(tmax),
           time = as.numeric(time)) %>% 
    mutate(time = time + tmax,
           time = log(time))
  
  fit <- gam_nb_regress(patdat)
  
  pred_data <- tibble(
    time = seq( min(patdat$time), log(90), length.out=1000 ),
    Method = "Discrete",
    delta = 1,
    vol = 1,
    disp_fct = 1,
    t_G = time+7,
    parentFraction = 1,
    bpr = 1) %>% 
    bind_rows(
      tibble(
        time = seq( min(patdat$time), log(10), length.out=100 ),
        Method = "Continuous",
        delta = 1,
        vol = 1,
        disp_fct = 1,
        t_G = time+7,
        parentFraction = 1,
        bpr = 1)) %>% 
    arrange(time)
  
  
    # Now predict
    fitted <- gratia::fitted_values(fit, data=pred_data, scale="link")
    
    pred_data <- pred_data %>% 
      left_join(fitted) %>% 
      # And add raw data
      ungroup() %>% 
      bind_rows(patdat) %>% 
      select(-patient) %>% 
      mutate(time = exp(time))
  
  return(pred_data)
  
  
}

preddata <- data %>% 
  group_by(patient) %>% 
  mutate(preds = map(patient, ~create_preds(data, .x)))

set.seed(12345)
measurements <- sample(unique(preddata$patient), 4, replace=F)




plotdata <- preddata %>% 
  ungroup() %>% 
  filter(patient %in% measurements) %>% 
  select(patient, preds) %>% 
  unnest(preds)

plotdata_discrete <- plotdata %>% 
  filter(Method=="Discrete")

plotdata_continuous <- plotdata %>% 
  filter(Method=="Continuous")

nb_indiv_fitplot <- ggplot(plotdata_continuous, aes(x=time, y=log(calc_aif))) +
  geom_point(size=1, shape=1, colour="grey") +
  geom_point(data=plotdata_discrete, colour="black", size=2) +
  #geom_point(data=plotdata_discrete, colour="red") +
  geom_line(data=plotdata_discrete %>% filter(!is.na(fitted)), 
            colour="red", aes(y=fitted)) +
  geom_line(data=plotdata_continuous %>% filter(!is.na(fitted)), 
            aes(y=fitted), linetype="dashed",colour="red") +
  facet_wrap(~patient, scales="free") +
  geom_ribbon(data=plotdata_discrete %>% filter(!is.na(fitted)), 
              aes(ymax=upper, ymin=lower), 
              alpha=0.3, fill="red" ) +
  theme_light() +
  labs(x="time (min)", y="log(AIF)") +
  NULL

nb_indiv_fitplot

ggsave(nb_indiv_fitplot, filename = "../Figures/nb_indiv_fitplot.png",
       width=6, height=5)
  
```



## 4. A model across all patients

Having resolved this issue of wiggliness, we wanted to improve our model further. Especially since our model is no longer shape-constrained, another way to incorporate more conservativism in the model is to model all the individuals at once in a hierarchical model.

Firstly, we found the mean of tmax for all patients and added back tmax for all patients' data.

```{r, message = FALSE, warning=FALSE}
table = table %>% mutate(
  Tmax = as.numeric(Tmax),
  Number_of_data = as.numeric(Number_of_data),
  sum = Tmax*Number_of_data
)

mean_tmax = sum(table$sum)/sum(table$Number_of_data)
```

```{r, message = FALSE, warning=FALSE}
total_data = data_frame()

for(i in 1:nrow(data)){
patient = data$patient[i] 
data_line = data$count_dsc[i] %>% as.data.frame()
p_data = cbind(patient, data_line)
total_data = rbind(total_data, p_data)
}

total_data = total_data %>% mutate(
  patient = as.factor(patient),
  time = log(time+mean_tmax)
  )

```

```{r, message = FALSE, warning=FALSE}
gam_total_nb_regress = function(data = data){
fit_res = gam(count ~ Method*patient+patient+s(time,k = 20)+s(time, patient, k=5, bs="fs"),
               offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+
                  (-log(calibration)) +
                  (-log(parentFraction))+log(bpr),
              family = nb(link = "log"),data = data,method="REML")
return(fit_res)
}

total_line = total_data %>% gam_total_nb_regress()
plot(total_line, main = "Regression for all patients")
qq.gam(total_line,rep=500)
 mtext(paste0("Q-Q plot for all measurements"), side = 3, line = -1, outer = TRUE)
```

GJM: For consistency, it would be nice to see plots like those you show above for the other models too. I think the QQ plot and individual smooths are so simple with one panel for everyone that it's worth keeping them, but I think it would be nice to show the individual predicted values to get an idea for whether this approach is underfitting the data.

```{r, message = FALSE, warning=FALSE}
create_preds_hgam <- function(data, patient_id) {
  
  patdat <- total_data %>% 
    filter(patient==patient_id) %>% 
    mutate(tmax = as.numeric(tmax),
           time = as.numeric(time)) %>% 
    mutate(time = time + tmax,
           time = log(time))
  
  pred_data <- tibble(
    time = seq( min(patdat$time), log(90), length.out=1000 ),
    Method = "Discrete",
    delta = 1,
    vol = 1,
    disp_fct = 1,
    t_G = time+7,
    parentFraction = 1,
    bpr = 1) %>% 
    bind_rows(
      tibble(
        time = seq( min(patdat$time), log(10), length.out=100 ),
        Method = "Continuous",
        delta = 1,
        vol = 1,
        disp_fct = 1,
        t_G = time+7,
        parentFraction = 1,
        bpr = 1)) %>% 
    arrange(time) %>% 
    mutate(patient = patient_id)
  
  
    # Now predict
    fitted <- gratia::fitted_values(total_line, 
                                    data=pred_data, scale="link")
    
    pred_data <- pred_data %>% 
      left_join(fitted) %>% 
      # And add raw data
      ungroup() %>% 
      bind_rows(patdat) %>% 
      select(-patient) %>% 
      mutate(time = exp(time))
  
  return(pred_data)
  
  
}



plot_hgam_preds <- function(patient, data, preds) {

  
  data_continuous <- data %>% 
    filter(Method=="Continuous")
  
  data_discrete <- data %>% 
    filter(Method=="Discrete")
  
  preds_continuous <- preds %>% 
    filter(Method=="Continuous") %>% 
    mutate(time = log(time))
  
  preds_discrete <- preds %>% 
    filter(Method=="Discrete") %>% 
    mutate(time = log(time))
  
  ggplot(data_continuous, aes(x=exp(time), y=log(calc_aif))) +
    geom_point(size=1, shape=1, colour="grey") +
    geom_point(data=data_discrete, colour="black", size=2) +
    #geom_point(data=plotdata_discrete, colour="red") +
    geom_line(data=preds_discrete, 
              colour="red", aes(y=fitted)) +
    geom_line(data=preds_continuous,
              aes(y=fitted), linetype="dashed",colour="red") +
    geom_ribbon(data=preds_discrete, 
                aes(ymax=upper, ymin=lower), 
                alpha=0.3, fill="red" ) +
    theme_light() +
    labs(x="time (min)", y="log(AIF)",
         subtitle = patient) +
    NULL
}

total_preddata <- total_data %>% 
  as_tibble() %>% 
  group_by(patient) %>% 
  nest() %>% 
  mutate(preds = map(patient, ~create_preds_hgam(total_data, .x)))

pmap(list(total_preddata$patient,
          total_preddata$data,
          total_preddata$preds),
     plot_hgam_preds)
```



## 5. Hierarchical extrapolation

### First and last samples

```{r, message = FALSE, warning=FALSE}
gam_total_nb_regress_missing = function(data = data, extrap_meas, total_mod){
  
  total_data <- data %>% 
    arrange(patient, time)
  
  extrap_id <- unique(data$patient)[extrap_meas]
  extrap_data <- data %>% filter(patient == extrap_id)
  
  disc_data <- extrap_data %>% 
    filter(Method=="Discrete")
  
  kept_data <- bind_rows(head(disc_data,1),
                         tail(disc_data,1))
  
  full_extrap_data <- extrap_data %>% 
    rename(calc_aif_hidden = calc_aif)
  
  meas_data <- data %>% anti_join(extrap_data) %>% 
    bind_rows(kept_data)
  
  pred_data = data_frame(
     patient = extrap_id,
     time = log(
       seq( 
         exp(extrap_data$time[1]), 
         90, 
         length.out=1000)),
     Method = "Discrete",
     delta = 1,
     vol = 1,
     disp_fct = 1,
     t_G = time+7,
     parentFraction = 1,
     bpr = 1
   )
  
  
  # Individual Model
  
  fit_res_indiv <- gam_nb_regress(extrap_data)
  
  pred_data$fit_indiv = predict(fit_res_indiv, newdata = pred_data)
  
  # Full Model
  
  pred_data$fit_interp = predict(total_mod, newdata = pred_data)
  
  # Extrap
  fit_res_extrap = gam(count ~ Method*patient+patient+s(time,k = 20)+s(time, patient, k=5, bs="fs"),
                 offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+
                    (-log(calibration)) +
                    (-log(parentFraction))+log(bpr),
                family = nb(link = "log"),data = meas_data,method="REML")
  
  pred_data$fit_extrap = predict(fit_res_extrap, newdata = pred_data)
  
  pred_data %>% 
    bind_rows(kept_data %>% select(patient, time, calc_aif)) %>% 
    bind_rows(full_extrap_data %>% select(patient, time, Method, calc_aif_hidden))
  
}

extrap_models_fl <- map(1:10, 
                     ~gam_total_nb_regress_missing(total_data, .x, total_line))

# for(i in 1:10) {
#   
#   pdat <- extrap_models[[i]]
#   
#   plot(pdat$time, pdat$fit_interp, type = "l")
#   lines(pdat$time, pdat$fit_extrap, type = "l", col="red")
#   points(pdat$time, log(pdat$calc_aif))
#   
# }

extrap_tibble_fl <- tibble(
  data = extrap_models_fl
) %>% 
  unnest()

extrap_fl <- ggplot(extrap_tibble_fl, aes(x=exp(time), y=log(calc_aif))) +
  geom_point(data = extrap_tibble_fl %>% filter(Method=="Discrete"),
             aes(y=log(calc_aif_hidden)), shape=1, colour="grey") +
  geom_point() +
  facet_wrap(~patient) +
  geom_line(data = extrap_tibble_fl %>% filter(!is.na(fit_extrap)),
            aes(y=fit_extrap), 
            colour="red") +
  geom_line(data = extrap_tibble_fl %>% filter(!is.na(fit_indiv)),
            aes(y=fit_indiv), colour="black", linetype="dotted") +
  theme_light() +
  labs(y="log(AIF)", x="time(min)") +
  NULL
  
extrap_fl
  
ggsave(extrap_fl, filename = "../Figures/Extrap_fl.png", width = 6, height=5)

# total_line = total_data %>% gam_total_nb_regress()
# plot(total_line, main = "Regression for all patients")
# qq.gam(total_line,rep=500)
#  mtext(paste0("Q-Q plot for all measurements"), side = 3, line = -1, outer = TRUE)
```

### Last samples

```{r, message = FALSE, warning=FALSE}
gam_total_nb_regress_missing = function(data = data, extrap_meas, total_mod){
  
  total_data <- data %>% 
    arrange(patient, time)
  
  extrap_id <- unique(data$patient)[extrap_meas]
  extrap_data <- data %>% filter(patient == extrap_id)
  
  disc_data <- extrap_data %>% 
    filter(Method=="Discrete")
  
  # kept_data <- bind_rows(head(disc_data,1),
  #                        tail(disc_data,1))
  
  kept_data <- tail(disc_data,2)
  
  full_extrap_data <- extrap_data %>% 
    rename(calc_aif_hidden = calc_aif)
  
  meas_data <- data %>% anti_join(extrap_data) %>% 
    bind_rows(kept_data)
  
  pred_data = data_frame(
     patient = extrap_id,
     time = log(
       seq( 
         exp(extrap_data$time[1]), 
         90, 
         length.out=1000)),
     Method = "Discrete",
     delta = 1,
     vol = 1,
     disp_fct = 1,
     t_G = time+7,
     parentFraction = 1,
     bpr = 1
   )
  
  
  # Individual Model
  
  fit_res_indiv <- gam_nb_regress(extrap_data)
  
  pred_data$fit_indiv = predict(fit_res_indiv, newdata = pred_data)
  
  # Full Model
  
  pred_data$fit_interp = predict(total_mod, newdata = pred_data)
  
  # Extrap
  fit_res_extrap = gam(count ~ Method*patient+patient+s(time,k = 20)+s(time, patient, k=5, bs="fs"),
                 offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+(-log(0.003))+(-0.0807*vol)
                 +(-log(parentFraction))+log(bpr)
                ,family = nb(link = "log"),data = meas_data,method="REML")
  
  pred_data$fit_extrap = predict(fit_res_extrap, newdata = pred_data)
  
  pred_data %>% 
    bind_rows(kept_data %>% select(patient, time, calc_aif)) %>% 
    bind_rows(full_extrap_data %>% select(patient, time, Method, calc_aif_hidden))
  
}

extrap_models_ll <- map(1:10, 
                     ~gam_total_nb_regress_missing(total_data, .x, total_line))

# for(i in 1:10) {
#   
#   pdat <- extrap_models[[i]]
#   
#   plot(pdat$time, pdat$fit_interp, type = "l")
#   lines(pdat$time, pdat$fit_extrap, type = "l", col="red")
#   points(pdat$time, log(pdat$calc_aif))
#   
# }

extrap_tibble_ll <- tibble(
  data = extrap_models_ll
) %>% 
  unnest()

extrap_ll <- ggplot(extrap_tibble_ll, aes(x=exp(time), y=log(calc_aif))) +
  geom_point(data = extrap_tibble_ll %>% filter(Method=="Discrete"),
             aes(y=log(calc_aif_hidden)), shape=1, colour="grey") +
  geom_point() +
  facet_wrap(~patient) +
  geom_line(data = extrap_tibble_ll %>% filter(!is.na(fit_extrap)),
            aes(y=fit_extrap), 
            colour="red") +
  geom_line(data = extrap_tibble_ll %>% filter(!is.na(fit_indiv)),
            aes(y=fit_indiv), colour="black", linetype="dotted") +
  theme_light() +
  labs(y="log(AIF)", x="time(min)") +
  NULL

extrap_ll
  
ggsave(extrap_ll, filename = "../Figures/Extrap_ll.png", width = 6, height=5)  

# total_line = total_data %>% gam_total_nb_regress()
# plot(total_line, main = "Regression for all patients")
# qq.gam(total_line,rep=500)
#  mtext(paste0("Q-Q plot for all measurements"), side = 3, line = -1, outer = TRUE)
```

### First samples

```{r, message = FALSE, warning=FALSE}
gam_total_nb_regress_missing = function(data = data, extrap_meas, total_mod){
  
  total_data <- data %>% 
    arrange(patient, time)
  
  extrap_id <- unique(data$patient)[extrap_meas]
  extrap_data <- data %>% filter(patient == extrap_id)
  
  disc_data <- extrap_data %>% 
    filter(Method=="Discrete")
  
  # kept_data <- bind_rows(head(disc_data,1),
  #                        tail(disc_data,1))
  
  kept_data <- head(disc_data,2)
  
  full_extrap_data <- extrap_data %>% 
    rename(calc_aif_hidden = calc_aif)
  
  meas_data <- data %>% anti_join(extrap_data) %>% 
    bind_rows(kept_data)
  
  pred_data = data_frame(
     patient = extrap_id,
     time = log(
       seq( 
         exp(extrap_data$time[1]), 
         90, 
         length.out=1000)),
     Method = "Discrete",
     delta = 1,
     vol = 1,
     disp_fct = 1,
     t_G = time+7,
     parentFraction = 1,
     bpr = 1
   )
  
  
  # Individual Model
  
  fit_res_indiv <- gam_nb_regress(extrap_data)
  
  pred_data$fit_indiv = predict(fit_res_indiv, newdata = pred_data)
  
  # Full Model
  
  pred_data$fit_interp = predict(total_mod, newdata = pred_data)
  
  # Extrap
  fit_res_extrap = gam(count ~ Method*patient+patient+s(time,k = 20)+s(time, patient, k=5, bs="fs"),
                 offset = log(delta)+log(vol)+log(disp_fct)+(-log(2)/20.364*t_G)+(-log(0.003))+(-0.0807*vol)
                 +(-log(parentFraction))+log(bpr)
                ,family = nb(link = "log"),data = meas_data,method="REML")
  
  pred_data$fit_extrap = predict(fit_res_extrap, newdata = pred_data)
  
  pred_data %>% 
    bind_rows(kept_data %>% select(patient, time, calc_aif)) %>% 
    bind_rows(full_extrap_data %>% select(patient, time, Method, calc_aif_hidden))
  
}

extrap_models_ff <- map(1:10, 
                     ~gam_total_nb_regress_missing(total_data, .x, total_line))

# for(i in 1:10) {
#   
#   pdat <- extrap_models[[i]]
#   
#   plot(pdat$time, pdat$fit_interp, type = "l")
#   lines(pdat$time, pdat$fit_extrap, type = "l", col="red")
#   points(pdat$time, log(pdat$calc_aif))
#   
# }

extrap_tibble_ff <- tibble(
  data = extrap_models_ff
) %>% 
  unnest()

extrap_ff <- ggplot(extrap_tibble_ff, aes(x=exp(time), y=log(calc_aif))) +
  geom_point(data = extrap_tibble_ff %>% filter(Method=="Discrete"),
             aes(y=log(calc_aif_hidden)), shape=1, colour="grey") +
  geom_point() +
  facet_wrap(~patient) +
  geom_line(data = extrap_tibble_ff %>% filter(!is.na(fit_extrap)),
            aes(y=fit_extrap), 
            colour="red") +
  geom_line(data = extrap_tibble_ff %>% filter(!is.na(fit_indiv)),
            aes(y=fit_indiv), colour="black", linetype="dotted") +
  theme_light() +
  labs(y="log(AIF)", x="time(min)") +
  NULL

extrap_ff
  
ggsave(extrap_ff, filename = "../Figures/Extrap_ff.png", width = 6, height=5)
  

# total_line = total_data %>% gam_total_nb_regress()
# plot(total_line, main = "Regression for all patients")
# qq.gam(total_line,rep=500)
#  mtext(paste0("Q-Q plot for all measurements"), side = 3, line = -1, outer = TRUE)
```

### Figure

```{r, fig.height=3, fig.width=6, message = FALSE, warning=FALSE}
set.seed(12345)

measurements <- sample(unique(total_data$patient), 4, replace=F)

extrap_combined <- bind_rows(
  extrap_tibble_fl %>% mutate(Samples = "First and last"),
  extrap_tibble_ff %>% mutate(Samples = "First two")
) %>% 
  filter(patient %in% measurements)

extrap_combined_plot <- ggplot(extrap_combined, aes(x=exp(time), y=log(calc_aif))) +
  geom_point(data = extrap_combined %>% filter(Method=="Discrete"),
             aes(y=log(calc_aif_hidden)), shape=1, colour="grey") +
  geom_point() +
  facet_grid( Samples ~ patient) +
  geom_line(data = extrap_combined %>% filter(!is.na(fit_extrap)),
            aes(y=fit_extrap), colour="red") +
  geom_line(data = extrap_combined %>% filter(!is.na(fit_extrap)),
            aes(y=fit_indiv), colour="black", linetype="dotted") +
  theme_light() +
  labs(y="log(AIF)", x="time (min)") +
  NULL


extrap_combined_plot

ggsave(extrap_combined_plot, filename = "../Figures/Extrap_plots.png", 
       width=6, height=3)
```

