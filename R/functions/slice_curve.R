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

