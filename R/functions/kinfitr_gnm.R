source("functions/gnm_wrap_decay.R")
# fit non-linear poisson model and return the model results
kinfitr_gnm <- function(t, # sampling time 
                        t_G = NULL, # time that put in the gamma counter
                        y.sum,  # count or concentration 
                        delta, # measure time in gamma counter
                        vol, # volume
                        decay_corrected = TRUE, # concentration decay corrected?
                        meta_corrected = TRUE, # concentration metabolism corrected?
                        pf = NULL, # fitted parent fractions. Later can be changed to a vector of parameters of the parent fraction function
                        bpr = NULL,
                        disp){
    # if input has been decay corrected
    if(decay_corrected == TRUE & !is.null(t_G)){
      # fit the model count ~ -1 + offset(log(delta)) + offset(log(vol)) + offset(-log(2)/20.364*t_G) + log_sum_exp(t) + offset(-log(0.003)*rep(1,length(t))) # calibration
      model_res = gnm_prop(t,t_G,y.sum,pf = rep(1,length(t)),bpr=rep(1,length(t)),disp,delta,vol)
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
      model_res = gnm_prop(t,t_G,y.sum,pf=pf,bpr=rep(1,length(t)),disp,delta,vol)

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
  