library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters,N) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I/N
    dI <-  beta * S * I/N - gamma * I/N
    dR <-                 gamma * I/N
    
    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1*10^6/(1*10^6), I = 100/(1*10^6), R = 0.0,N=1*10^6)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)
## Time frame
times      <- seq(0, 20, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)
plot(out$I,type='l')
