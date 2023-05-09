from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Define initial population sizes and parameters
N = 10000   # total population
E0, I0, R0 = 10, 1, 0   # initial numbers of exposed, infected, and recovered individuals
S0 = N - E0 - I0 - R0   # initial number of susceptible individuals
beta = 0.4   # average number of people an infected person infects per day
gamma = 0.1   # fraction of infected individuals that recover per day
sigma = 0.2   # fraction of susceptible individuals that become exposed per day
t = np.linspace(0, 100, 100)   # time interval in days

# Define SEIR model function


def SEIR_model(y, t, N, beta, gamma, sigma):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt


# Solve SEIR model
y0 = S0, E0, I0, R0
sol = odeint(SEIR_model, y0, t, args=(N, beta, gamma, sigma))
S, E, I, R = sol.T

# Plot results
plt.plot(t, S, 'b', label='Susceptible')
plt.plot(t, E, 'y', label='Exposed')
plt.plot(t, I, 'r', label='Infected')
plt.plot(t, R, 'g', label='Recovered')
plt.legend(loc='best')
plt.xlabel('Time (days)')
plt.ylabel('Population')
plt.title('SEIR Model')
plt.show()
