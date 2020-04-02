library(rjags)
library(R2jags)
library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = .95, gamma = 1/3.5)
## Time frame
times      <- seq(0, 35, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
plot(out$I)

model <- "
model {
   
    beta ~ dnorm(.95,100000000)
    gamma ~ dnorm(1/3.5,100000000)
    I_prop_init ~ dbeta(100000*1e-6,100000*(1-1e-6))
    S[1] <- round(n_sample*(1-I_prop_init))
    I[1] <- round(n_sample*I_prop_init)
    R[1] <- 0
    delta_t ~ dbeta(10*.5,10*(1-.5))
    for (i in 2:n_obs){
    
        delta_n_si[i] ~ dbinom(1-exp(-beta*I[i-1]/n_sample*delta_t),S[i-1])
        delta_n_ir[i] ~ dbinom(1-exp(-gamma*delta_t),I[i-1])
        S[i] <- S[i-1] - delta_n_si[i]
        I[i] <- I[i-1] + delta_n_si[i] - delta_n_ir[i]
        R[i] <- R[i-1] + delta_n_ir[i]
    }

    for (i in 1:n_obs){
      y[i] ~ dnorm(I[i],1000000)
    }
    #for (j in 14:d_length){
    #    d[j] ~ dnorm(I[j-13]*.0134 ,1000)
   # }

}"
state_testing_data <- read.csv("/Users/gcgibson/Downloads/states-daily.csv")

state_testing_data$date_formatted <- as.Date(state_testing_data$dateChecked)

state_testing_data$date_formatted <- unlist(lapply(state_testing_data$date,function(x){
  year <- substr(x,1,4)
  month <- substr(x,5,6)
  day <- substr(x,7,8)
  return (paste0(year,"/",month,"/",day))
}))
state_testing_data$date_formatted <- as.Date(state_testing_data$date_formatted)
state_testing_data$week <- lubridate::week(state_testing_data$date_formatted)

state_testing_data_by_week <- state_testing_data %>% group_by(week,state) %>% summarise(week_positive = sum(positive,na.rm=T))
state_testing_data <- state_testing_data %>% arrange(state_testing_data$date_formatted)
state_testing_data_ny <- state_testing_data[state_testing_data$state=="NY",]
# For stan model we need the following variables:
state_testing_data_ny$death[is.na(state_testing_data_ny$death)] <- 0

# subset to last 20 rows
nrow_df <- nrow(state_testing_data_ny)
state_testing_data_ny_subset <- state_testing_data_ny[tail(1:nrow_df,20),]


stan_d = list(n_obs = length(state_testing_data_ny_subset$positive)+7*8,
              n_params = 3,
              n_difeq = 4,
              n_sample = 10^6,
              y = c(state_testing_data_ny_subset$positive,rep(NA,7*8)),
              #y=rep(NA,length(state_testing_data_ny_subset$positive)+7*4),
              d = state_testing_data_ny_subset$death,
              d_length = length(state_testing_data_ny_subset$death))


jgs <- jags.model(file = textConnection(model), data = stan_d, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('I','beta','gamma'), 3000, 3)
plot(rowMeans(out_jags$I))
lines(state_testing_data_ny_subset$positive,col='red')
#hist(out_jags$beta)
#hist(out_jags$gamma)
