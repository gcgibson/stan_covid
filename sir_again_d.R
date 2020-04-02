library(deSolve)
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

cat(
  '
  functions {
  
  // This largely follows the deSolve package, but also includes the x_r and x_i variables.
  // These variables are used in the background.
  
  real[] SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  
  real dydt[4];
  
  dydt[1] = - params[1] * y[1] * y[2];
  dydt[2] = params[1] * y[1] * y[2] - params[3] * y[2];
  dydt[3] =  params[3] * y[2] - params[2] * y[3];
  dydt[4] = params[2] * y[3];

  return dydt;
  }
  
  }
  
  data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_sample; // Number of hosts sampled at each time point.
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  
  real y[n_obs]; // The binomially distributed data
  int d[n_obs]; // The binomially distributed data
  
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  }
  
  parameters {
  real<lower = 0> params[n_params]; // Model parameters
  real<lower = 0, upper = 1> S0; // Initial fraction of hosts susceptible
  real y_latent[n_obs];
  }
  
  transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real y0[n_difeq]; // Initial conditions for both S and I
  

  y0[1] = S0;
  y0[2] = (1 - S0)/2;
  y0[3] = (1-S0)/2;
  y0[4] = 0;

  y_hat = integrate_ode_rk45(SI, y0, t0, ts, params, x_r, x_i);
  
  }
  
  model {
  // parameter estimates come from
  // https://www.medrxiv.org/content/10.1101/2020.03.21.20040303v2.full.pdf+html
  params[1] ~ normal(.95, .001); //prior estimate for beta from 
  params[2] ~ normal(1/3.5, .001); //constrained to be positive
  params[3] ~ normal(1/3.65, .001); //constrained to be positive

  S0 ~ normal(0.5, 0.5); //constrained to be 0-1.
  for (i in 1:n_obs){
    y_latent[i] ~ normal(n_sample*y_hat[i, 2],.0001); //y_hat[,2] are the fractions infected from the ODE solver
  }
  y ~ normal(y_latent,1);
   # for (d_i in 14:n_obs){
   #   d[d_i] ~ normal(y[n_obs-13]*.0134,1); //y_hat[,2] are the fractions infected from the ODE solver
   # }
  }
  
  generated quantities {
  // Generate predicted data over the whole time series:
  real fake_I[n_fake, n_difeq];
  
  fake_I = integrate_ode_rk45(SI, y0, t0, fake_ts, params, x_r, x_i);
  
  }
  
  ', 
  file = "SI_fit.stan", sep="", fill=T)

# FORMAT Data
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
stan_d = list(n_obs = length(state_testing_data_ny$positive),
              n_params = 3,
              n_difeq = 4,
              n_sample = 10^6,
              n_fake = length(1:length(state_testing_data_ny$positive))+7*4,
              y = state_testing_data_ny$positive,
              d = state_testing_data_ny$death,
              t0 = 0,
              ts = 1:length(state_testing_data_ny$positive),
              fake_ts = c(1:(length(state_testing_data_ny$positive)+7*4)))

# Which parameters to monitor in the model:
params_monitor = c("y_hat", "y0", "params", "fake_I")

# Test / debug the model:
test = stan("SI_fit.stan",
            data = stan_d,
            pars = params_monitor,
            chains = 1, iter = 10)

# Fit and sample from the posterior
mod = stan(fit = test,
           data = stan_d,
           pars = params_monitor,
           chains = 1,
           warmup = 500,
           iter = 1500)

# You should do some MCMC diagnostics, including:
#traceplot(mod, pars="lp__")
#traceplot(mod, pars=c("params", "y0"))
#summary(mod)$summary[,"Rhat"]

# These all check out for my model, so I'll move on.

# Extract the posterior samples to a structured list:
posts <- extract(mod)
plot(colMeans(posts$fake_I[,,3])*10^6)
library(ggplot2)
col_num <- ncol(posts$fake_I)
df_for_plot <- data.frame(x=rep(1:col_num,length.out=3000*col_num),y=c(t(posts$fake_I[,,3])))
ggplot(df_for_plot,aes(x=x,y=y*10^6)) + geom_point(alpha=.1) + theme_bw() + geom_line(data=state_testing_data_ny,aes(x=1:length(positive),y=positive,col='truth'))+
  geom_line(data=state_testing_data_ny,aes(x=1:length(death),y=death,col='truth'))

