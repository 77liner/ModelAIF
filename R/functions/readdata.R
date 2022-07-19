library(tidyverse)
# read data from 10 subjects
readdata <- function(data_type = NULL){
    filenames = list.files(path = "./data", pattern="*.csv") ##name of files
    data = tibble(
      patient = substr(filenames,1,5), #the first 5 character of filenames
      count = map(filenames,~read.csv(paste0("./data/", .)) %>% select(-X))
    )
  return(data)
}