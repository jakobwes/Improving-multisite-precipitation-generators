dat <- read_csv("data/2_data_interim/interim_data.csv")

get_resolution <- function(precip){
  precip_sorted <- sort(precip)
  differences <- precip_sorted[2:length(precip_sorted)] - precip_sorted[1:(length(precip_sorted)-1)]
  
  return(
    min(differences[differences != 0])
  )
}

# Plot yearly resolution (minimum difference in values)
dat %>% mutate(year = year(date)) %>% group_by(site, year) %>% 
  summarise(resolution = get_resolution(precip)) %>%
  ggplot(aes(year, resolution)) + geom_point() + ylim(0, 1) + facet_wrap(~site)

# We see a changepoint at around 1970. Investigate that closer for one site
dat %>% filter(site == "S003") %>% mutate(year = year(date)) %>% group_by(year) %>% 
  summarise(resolution = get_resolution(precip)) %>%
  filter(year < 1980) %>%
  ggplot(aes(year, resolution)) + geom_point() + ylim(0, 1) 

# Some years have an even higher resolution than 0.1. There seems to be no systematic structure to it, but investigate how many values this concerns
dat %>% mutate(diff = abs(round(precip, digits = 1) - precip)) %>% filter(diff > 0.0001) %>% pull(diff)
