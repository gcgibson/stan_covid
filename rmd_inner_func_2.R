run_model_for_location <- function(location){
  model <- "
  model {
  
  beta ~ dnorm(.95,1000)
  gamma ~ dnorm(1/3.5,1000)
  #sigma ~ dnorm(1/3.67,1000)
  #beta <- 10
  #gamma <- .5
  
  S[1] <- n_sample - 1
  I[1] <- 1
  R[1] <- 0
  t[1] <- 1
  for (i in 2:1000){
      w1[i] <- beta*S[i-1]*I[i-1]/1e3
      w2[i] <- gamma*I[i-1]
      W[i] <- w1[i] + w2[i]
      U[i] ~ dunif(0,1)
      dt[i] <- -log(U[i]) /W[i]
      t[i] <- t[i-1] + dt[i]

      U2[i] ~ dunif(0,1)
      s1[i] <- ifelse(U2[i] < (w1[i] / W[i]),-1,1)
      s2[i] <- ifelse(U2[i] < (w1[i] / W[i]),0,1)
        

      S[i] <- S[i-1] + s1[i]*(1-s2[i])  
      I[i] <- max(I[i-1] - s1[i],1)
  }
  
  
  for (i in 1:n_obs){
      y[i] ~ dnorm(I[i*10],10)
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
                y1 = state_testing_data_ny_subset$positive[1],
                #y=rep(NA,length(state_testing_data_ny_subset$positive)+h),
                d = state_testing_data_ny_subset$death,
                d_length = length(state_testing_data_ny_subset$death))
  
  library(R2jags)
  library(rjags)
  jgs <- jags.model(file = textConnection(model), data = stan_d, n.adapt = 1000)
  update(jgs, 1000)
  out_jags <- jags.samples(jgs, c('I','beta','gamma','y','t'), 3000, 3)
  plot(rowMeans(out_jags$I))
  rowMeans(out_jags$t)
  #lines(state_testing_data_ny_subset$positive,col='red')
  library(ggplot2)
  data_for_plot <- data.frame(x=rep(1:dim(out_jags$y)[1],dim(out_jags$y)[2]),
                              y=c(out_jags$y),group=rep(1:dim(out_jags$y)[2],each=dim(out_jags$y)[1]))
  
  forecast_plot <- ggplot(data_for_plot,aes(x=x,y=y,group=group)) + geom_line(alpha=.1) + theme_bw() + geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$positive),aes(x=x,y=y,group=1,col='observed')) 
  
  nowcast_plot <- ggplot(data_for_plot,aes(x=x,y=y,group=group)) + geom_line(alpha=.1) + theme_bw() + geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$positive),aes(x=x,y=y,group=1,col='observed')) +
    geom_line(data=data.frame(x=1:length(state_testing_data_ny_subset$positive),y=state_testing_data_ny_subset$death),aes(x=x,y=y,group=1,col='death observed'))
  
  beta_hist <- ggplot(data=data.frame(x=as.vector(out_jags$beta))) + geom_histogram(aes(x=x)) + xlab("beta")
  
  gamma_hist <- ggplot(data=data.frame(x=as.vector(out_jags$gamma))) + geom_histogram(aes(x=x)) + xlab("gamma")
  p_grid <- plot_grid(beta_hist,gamma_hist,forecast_plot,ncol=1,rel_heights = c(1,1,2))
  return (p_grid)
}

run_model_for_location("NJ")
#print (run_model_for_location("OH"))
#print (run_model_for_location("MA"))
#print (run_model_for_location("NY"))

