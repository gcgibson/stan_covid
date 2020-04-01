library(FluSight)
library(cdcfluutils)
state_testing_data <- read.csv("/Users/gcgibson/Downloads/states-daily.csv")
counties_testing_data <- read.csv("/Users/gcgibson/Downloads/counties.csv")

library(ggplot2)
state_testing_data$date_formatted <- as.Date(state_testing_data$dateChecked)

state_testing_data$date_formatted <- unlist(lapply(state_testing_data$date,function(x){
  year <- substr(x,1,4)
  month <- substr(x,5,6)
  day <- substr(x,7,8)
  return (paste0(year,"/",month,"/",day))
}))
state_testing_data$date_formatted <- as.Date(state_testing_data$date_formatted)

library(R2jags)
library(rjags)

model <- "
model {
  #
  positive_latent[1] ~ dnorm(1000,.000001)
  for (i in 2:N){
    positive_latent[i] ~ dnorm(positive_latent[i-1],.000001)
  }
  
  # integer valued delay model based on previous studies
  onset_to_death ~ dbinom(.75,21)
  cfr ~ dbeta(1000*.013,1000*(1-.013))
  for (i in 28:N){
    deaths[i-27] ~ dbinom(cfr,round(positive_latent[i-min(28,onset_to_death)]))
  }
}"


state_testing_data_d_pos <- state_testing_data#[state_testing_data$death >0,]
state_testing_data_d_pos <- state_testing_data_d_pos %>% dplyr::arrange(date_formatted)
positive_cases <- state_testing_data_d_pos[state_testing_data_d_pos$state == "NY",]$positive
deaths <-state_testing_data_d_pos[state_testing_data_d_pos$state == "NY",]$death

dat <- list(positive = positive_cases,deaths=deaths,N=length(deaths)+27)


## add to dataframe
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('positive_latent'), 3000, 3)

rowMeans(out_jags$positive_latent)[28:length(rowMeans(out_jags$positive_latent))]
deaths
