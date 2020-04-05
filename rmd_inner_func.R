run_model_for_location <- function(location){
  model <- "
  model {
  
  beta ~ dnorm(.95,1000)
  gamma ~ dnorm(1/3.5,1000)
  sigma ~ dnorm(1/3.67,1000)
  
  I_prop_init ~ dbeta(100*1e-6,100*(1-1e-6))
  E_prop_init ~ dbeta(100*1e-6,100*(1-1e-6))
  
  S[1] <- round(n_sample*(1-I_prop_init-E_prop_init))
  E[1] <- round(n_sample*E_prop_init)
  I[1] <- round(n_sample*I_prop_init)
  R[1] <- 0
  delta_t ~ dbeta(100*.1,100*(1-.1))
  
  var_transition <- 10
  for (i in 2:n_obs){
      p_se[i] <- max(1-exp(-beta*I[i-1]/n_sample*delta_t),1e-10)
      p_ei[i] <- max(1-exp(-sigma*delta_t),1e-10)
      p_ir[i] <- max(1-exp(-gamma*delta_t),1e-10)
    
      delta_n_se[i] ~ dbeta(var_transition*p_se[i],var_transition*(1-p_se[i]))
      delta_n_ei[i] ~ dbeta(var_transition*p_ei[i],var_transition*(1-p_ei[i]))
      delta_n_ir[i] ~ dbeta(var_transition*p_ir[i],var_transition*(1-p_ir[i]))
      
      S[i] <- S[i-1] - delta_n_se[i]*S[i-1]
      E[i] <- E[i-1] + delta_n_se[i]*S[i-1] -delta_n_ei[i]*E[i-1]
      I[i] <- I[i-1] + delta_n_ei[i]*E[i-1] - delta_n_ir[i]*I[i-1]
      R[i] <- R[i-1] + delta_n_ir[i]*I[i-1]
  }
  
  for (i in 1:n_obs){
    y[i] ~ dnorm(I[i],10)
  }
  
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
  
  
  state_testing_data <- state_testing_data %>% arrange(state_testing_data$date_formatted)
  state_testing_data_ny <- state_testing_data[state_testing_data$state==location,]
  
  
  # subset to last 20 rows
  nrow_df <- nrow(state_testing_data_ny)
  start_idx <- which(state_testing_data_ny$positive > 2)[1]
  state_testing_data_ny_subset <- state_testing_data_ny[start_idx:nrow_df,]
  
  
  h <- 7*8
  stan_d = list(n_obs = length(state_testing_data_ny_subset$positive)+h,
                n_sample = 10^6,
                y = c(state_testing_data_ny_subset$positive,rep(NA,h)),
                #y=rep(NA,length(state_testing_data_ny_subset$positive)+h),
                d = state_testing_data_ny_subset$death,
                d_length = length(state_testing_data_ny_subset$death))
  
  library(R2jags)
  library(rjags)
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
  
  beta_hist <- ggplot(data=data.frame(x=as.vector(out_jags$beta))) + geom_histogram(aes(x=x)) + xlab("beta")
  
  gamma_hist <- ggplot(data=data.frame(x=as.vector(out_jags$gamma))) + geom_histogram(aes(x=x)) + xlab("gamma")
  p_grid <- plot_grid(beta_hist,gamma_hist,forecast_plot,ncol=1,rel_heights = c(1,1,2))
  return (p_grid)
}

print (run_model_for_location("NJ"))
print (run_model_for_location("OH"))
print (run_model_for_location("MA"))
print (run_model_for_location("NY"))

