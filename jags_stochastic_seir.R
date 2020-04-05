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
   
    beta ~ dnorm(.95,10000)
    gamma ~ dnorm(1/3.5,10000)
    sigma ~ dnorm(1/3.67,10000)

    #E_prop_init ~ dbeta(100000*1e-5,100000*(1-1e-5))
    I_prop_init ~ dbeta(100000*1e-5,100000*(1-1e-5))

    S[1] <- round(n_sample*(1-I_prop_init))
    E[1] <- 0
    I[1] <- round(n_sample*I_prop_init)
    R[1] <- 0
    delta_t ~ dbeta(100*.1,100*(1-.1))
    for (i in 2:n_obs){
    
        #delta_n_se[i] ~ dbinom(1-exp(-beta*I[i-1]/n_sample*delta_t),S[i-1])
        #delta_n_ei[i] ~ dbinom(1-exp(-sigma*delta_t),E[i-1])
        #delta_n_ir[i] ~ dbinom(1-exp(-gamma*delta_t),I[i-1])
        
        p_se[i] <- 1-exp(-beta*I[i-1]/n_sample*delta_t)
        p_ei[i] <- 1-exp(-sigma*delta_t)
        p_ir[i] <- 1-exp(-gamma*delta_t)
        delta_n_se[i] ~ dbeta(10*p_se[i],10*(1-p_se[i]))
        delta_n_ei[i] ~ dbeta(10*p_ei[i],10*(1-p_ei[i]))
        delta_n_ir[i] ~ dbeta(10*p_ir[i],10*(1-p_ir[i]))
        
        S[i] <- S[i-1] - delta_n_se[i]*S[i-1]
        E[i] <- E[i-1] + delta_n_se[i]*S[i-1] -delta_n_ei[i]*E[i-1]
        I[i] <- I[i-1] + delta_n_ei[i]*E[i-1] - delta_n_ir[i]*I[i-1]
        R[i] <- R[i-1] + delta_n_ir[i]*I[i-1]
    }

    for (i in 1:n_obs){
      y[i] ~ dnorm(I[i],1)
    }

#   death_dist ~ dcat(mPriorProb[])
  # for (prior_cat_probs in 1:5){
  #    mPriorProb[prior_cat_probs] <- 0
  # }
  
  # mPriorProb[6] <- .05
   #mPriorProb[7] <- .05
  # mPriorProb[8] <- .1
  # mPriorProb[9] <- .1
  # mPriorProb[10] <- .1
  # mPriorProb[11] <- .2
  # mPriorProb[12] <- .2
  # mPriorProb[13] <- .2
  # mPriorProb[14] <- .2
  # mPriorProb[15] <- .2
  # mPriorProb[16] <- .2
  # mPriorProb[17] <- .2
  # mPriorProb[18] <- .2
  # mPriorProb[19] <- .2
  # mPriorProb[20] <- .2
  # mPriorProb[21] <- .2
  # mPriorProb[22] <- .2
  # mPriorProb[23] <- .1
  # mPriorProb[24] <- .1



   #for (j in 14:d_length){
    #    idx[j] <- max(j-death_dist,1)
   #     d[j] ~ dnorm(I[idx[j]]*.0134 ,1000)
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
h <- 7*4

stan_d = list(n_obs = length(state_testing_data_ny_subset$positive)+h,
              n_sample = 10^6,
              y = c(state_testing_data_ny_subset$positive,rep(NA,h)),
              #y=rep(NA,length(state_testing_data_ny_subset$positive)+7*4),
              d = state_testing_data_ny_subset$death,
              d_length = length(state_testing_data_ny_subset$death))


jgs <- jags.model(file = textConnection(model), data = stan_d, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('I','beta','gamma','y'), 3000, 3)
#plot(rowMeans(out_jags$I))
#lines(state_testing_data_ny_subset$positive,col='red')
library(ggplot2)
data_for_plot <- data.frame(x=rep(1:dim(out_jags$I)[1],dim(out_jags$I)[2]),
                            y=c(out_jags$I),group=rep(1:dim(out_jags$I)[2],each=dim(out_jags$I)[1]))

forecast_plot <- ggplot(data_for_plot,aes(x=x,y=y,group=group)) + geom_line(alpha=.1) + theme_bw() + geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$positive),aes(x=x,y=y,group=1,col='observed'))
nowcast_plot <- ggplot(data_for_plot,aes(x=x,y=y,group=group)) + geom_line(alpha=.1) + theme_bw() + geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$positive),aes(x=x,y=y,group=1,col='observed')) +
  geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$death),aes(x=x,y=y,group=1,col='death observed'))


hist(out_jags$beta)
hist(out_jags$gamma)
