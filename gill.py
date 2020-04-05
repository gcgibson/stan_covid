import pandas as pd
import sys
import math
import random
import numpy as np
import scipy
import scipy.integrate as spi



def run_sim_location(loc,pop):
 
  
  
  def logit(p):
      return (math.log(p/(1-p)))
  
  def inv_logit(x):
      return (1/(1+math.exp(-x)))
  
  
  def forward_simulate_deterministic(delta_t,beta,gamma,sigma,E_init,I_init):
      def diff_eqs(INP,t):  
          Y=np.zeros((4))
          V = INP    
          Y[0] = - beta*(.93**t) * V[0] * V[2]
          Y[1] = beta*(.93**t) * V[0] * V[2] - sigma * V[1]
          Y[2] = sigma*V[1] - gamma*V[2]
          Y[3] = gamma * V[2]
          return Y   # For odeint
      INPUT = (1-E_init-I_init, E_init, I_init, 0.0)
      t_start = 0.0; t_end = len(confirmed_ny)+7*12; t_inc = delta_t
      t_range = np.arange(t_start, t_end+t_inc, t_inc)
      RES = spi.odeint(diff_eqs,INPUT,t_range)
      return (RES[:,2])

 
  
  
  def forward_simulate_stochastic(delta_t,beta,gamma):
      I_init = 10
      n_I = [I_init]
      n_S = [N - I_init]
      n_R = [0]
      for t in range(1,len(confirmed_ny)+40):
          if n_I == 0:
              break
          p_se = max(1-np.exp(-delta_t*beta*n_I[t-1]/N),1e-10)
          p_ei = max(1-np.exp(-delta_t*gamma),1e-10)
          logit_pse = logit(p_se)
          logit_pei = logit(p_ei)
          p_se = inv_logit(logit_pse +np.random.normal(0,.0001))
          p_ei = inv_logit(logit_pei +np.random.normal(0,.0001))
          n_S.append(n_S[t-1] - p_se*n_S[t-1])
          n_I.append(n_I[t-1] + p_se*n_S[t-1] - p_ei*n_I[t-1])
      return (n_I)
  
  
  
  def load_us():
      df = pd.read_csv('https://covidtracking.com/api/states/daily.csv')
      df.date = pd.to_datetime(df.date, format='%Y%m%d')
      df = df.set_index('date')
      df = df.drop(columns=['dateChecked'])
      df['confirmed'] = df['positive'] # useful alias to match other data
      df = df.pivot(columns='state')
      df.columns = df.columns.swaplevel()
      df.sort_index(axis=1, inplace=True)
      return df
  
  covid = load_us() 
  confirmed_ny = np.array(covid[loc]['confirmed'])[-10:]
  #print (covid['NY'])
  #sys.exit()
 
  N =pop
  
  
  I_init = 1e-5   
  E_init = 5*I_init
  
  
  # Initialize results list
  
  delta_t_vec = []
  beta_vec = []
  gamma_vec = []
  mse_vec = []
  for nsim in range(1000):
  # Main loop
      delta_t = 1
      var = 1
      beta = np.random.normal(.95,var)
      gamma = 1/3.5#np.random.normal(1/3.5,var)
      sigma = 1/3.7#np.random.normal(1/3.7,var)
  
      n_I = forward_simulate_deterministic(delta_t,beta,gamma,sigma,E_init,I_init)*N
  
      mse_vec.append(np.sum((n_I[0:len(confirmed_ny)]-confirmed_ny)**2)/len(confirmed_ny))
      #if (mse_vec[nsim] < mse_vec[nsim-1]): # or it doesn't draw of 
      delta_t_vec.append(delta_t)
      beta_vec.append(beta)
      gamma_vec.append(gamma)
  
  
  import matplotlib.pyplot as plt
  

  beta_optimal = beta_vec[mse_vec.index(min(mse_vec))]
  
  n_I_simulate_matrix = []
  for nsim_param in range(len(beta_vec)):
      n_I = forward_simulate_deterministic(1,beta_optimal,gamma_vec[nsim_param],sigma,E_init,I_init)
      n_I = N*np.random.beta(500*n_I,500*(1-n_I))
      n_I_simulate_matrix.append(n_I)
  
  
  plot= False
  if plot == True:
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.scatter(np.repeat(range(len(n_I)),len(beta_vec)),np.array(np.transpose(n_I_simulate_matrix)).flatten(), s=10, c='r', marker="o", label='predicted',alpha=.1)
      ax1.scatter(range(len(confirmed_ny)), confirmed_ny, s=10, c='b', marker="s", label='confirmed')
      plt.legend(frameon=False, loc="upper left")
      plt.tick_params(axis="both", which="major",labelsize=6,direction="in")
      #plt.style.use("fivethirtyeight");
      #plt.ylim(0,N);
      #plt.yscale("log")
      plt.show()
  
  return (n_I_simulate_matrix)


