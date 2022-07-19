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