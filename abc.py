import numpy as np
import numpy.matlib as nm
import scipy
import scipy.integrate as spi
from scipy.stats import norm,gamma
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize, rosen, rosen_der


def run_sim_location(loc,pop):

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


    def simulate(beta,gamma,I_init,N):

       
        sigma = 1/3.7
        def deriv(y, t, N, beta, gamma,sigma):
            S, E,I, R = y
            dSdt = -beta*.985**t * S * I / N
            dEdt = beta*.985**t * S * I / N - sigma * E
            dIdt = sigma*E - gamma*I
            dRdt = gamma * I

            return dSdt,dEdt, dIdt, dRdt

        # Initial conditions vector
        y0 = 1-4*I_init,3*I_init, I_init, 0
        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, y0, range(N), args=(1, beta, gamma,sigma))
        S, E,I, R = ret.T
        return (I)
    def run_abc(data,N):
        betas = np.random.normal(.95,.0001,size=1000)
        gammas = np.random.normal(1/3.5,.0001,size=1000)
        I_init = np.random.normal(.0001,.0001,size=1000)
        summary_stats = []
        for param_idx in range(1000):
            sims = simulate(betas[param_idx],gammas[param_idx],I_init[param_idx], N)
            if len(sims) == len(data):
                summary_stats.append(np.sum((sims-data)**2))
            else:
                summary_stats.append(1e20)
        return betas,gammas,I_init,summary_stats


    def log_prob(theta,data,N):
            ret_likelihood = (-1*(np.sum((simulate(theta,N)-data)**2)))
            ret_proir = -.5 *((theta-1)/1)**2
            ret_prob = ret_likelihood + ret_proir
            return ret_prob
    def dlog_prob(theta):
            #print(theta)
            ret_likelihood = (-1*(np.sum((simulate(theta,80)-sim))))
            ret_proir = -((theta-1)/1)
            ret_prob = ret_likelihood + ret_proir
            return ret_prob
    def dlog_prob_total(theta):
            total_log_prob = []
            for theta_l in theta:
                total_log_prob.append(dlog_prob(theta_l))
            return (np.array(total_log_prob).reshape((-1,1)))

    def kde_expansion(params,summary_stats,tol):
        possible_betas = []
        for i in range(len(summary_stats)):
            if (summary_stats[i]) < tol:
                possible_betas.append(params[i])

        kde_expansion = []
        for i in possible_betas:
            kde_expansion.append(np.random.normal(i,.0001,size=10000))
        return (np.array(kde_expansion).flatten())

    
    def mle(param):
        beta,gamma,i_init = param
        t_sol = simulate(beta,gamma,i_init, len(sim))
        return np.linalg.norm(t_sol-sim);


    import emcee

    covid = load_us() 
    confirmed_ny = np.array(covid[loc]['confirmed'])
    confirmed_ny_start_idx =np.nonzero(confirmed_ny)[0][0]

    confirmed_ny = confirmed_ny[confirmed_ny_start_idx:]
    confirmed_ny = np.array([x for x in confirmed_ny if str(x) != 'nan'])

    N_local = len(confirmed_ny)
    sim = confirmed_ny/(pop)
    res = minimize(mle, np.array([.95,1/3.5,.0001]),method= 'SLSQP')
    print (res)
    print (sim)
    #betas,gammas,I_init,summary_stats= run_abc(sim,len(confirmed_ny))
    #tol = .1
    #betas_expanded = kde_expansion(betas,summary_stats,tol)
    #gammas_expanded = kde_expansion(gammas,summary_stats,tol)
    #init_expanded = kde_expansion(I_init,summary_stats,tol)
    betas_expanded = np.random.normal(loc=res.x[0],scale=.00001,size=1000)
    gammas_expanded = np.random.normal(loc=res.x[1],scale=.00001,size=1000)
    init_expanded = np.random.normal(loc=res.x[2],scale=.00001,size=1000)
   
    n_I_simulate_matrix = []

    for i in range(1000):
        sim_res = simulate(betas_expanded[i],gammas_expanded[i],init_expanded[i],100)
        plt.plot(sim_res,color='blue',alpha=.1)
        n_I_simulate_matrix.append(sim_res)
    #plt.plot(sim,color='red')
    ##plt.ylim(0,.2)
    #plt.show()
    return (n_I_simulate_matrix)
