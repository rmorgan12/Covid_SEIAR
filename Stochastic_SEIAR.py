import numpy as np
import matplotlib.pyplot as plt


# Define the SEIR model parameters
total_pop = 220000000
m = .1/(67.7*365)  # natural mortality
a = m*total_pop  # recruitement Rate
b = .9549  # transmission rate
n = .9635  # transmissibility multiple of asymptomatic infection
s = .8961
w = .9692
m0 = .0126  # death rate during E(t)
m1 = .001  # death rate during I(t)
m2 = .0069  # death rate druing A(t)
y1 = .1456  # recovery rate for symp
y2 = .8666  # recovery rate during for asymp


# Define the initial values for the model variables
S0 = 219698286
E0 = 300000
I0 = 1714
A0 = 0
R0 = 0
N0 = total_pop

# Define the duration of the simulation (in days)
t_max = 365
dt = 1/365
sd_list = [0, .1, .2, .3]

# Define the time vector
t = np.arange(0, t_max, dt)
print(len(t))
# Define the arrays to store the model variables
S = np.zeros(len(t))
E = np.zeros(len(t))
I = np.zeros(len(t))
A = np.zeros(len(t))
R = np.zeros(len(t))
N = np.zeros(len(t))

# Initialize the arrays with the initial values
S[0] = S0
E[0] = E0
I[0] = I0
A[0] = A0
R[0] = R0
N[0] = N0
paths = []
np.random.seed(8)
# Implement Euler's method to simulate the SEIAR model
for sd in sd_list:
    S = np.zeros(len(t))
    E = np.zeros(len(t))
    I = np.zeros(len(t))
    A = np.zeros(len(t))
    R = np.zeros(len(t))
    N = np.zeros(len(t))
    S[0] = S0
    E[0] = E0
    I[0] = I0
    A[0] = A0
    R[0] = R0
    N[0] = N0

    for i in range(1, len(t)):
        dSdt = (a - b*(I[i-1] + n*A[i-1])*(S[i-1]/N[i-1]) - m *
                S[i-1])*dt + sd*np.sqrt(dt)*np.random.normal(0, 1)*S[i-1]
        dEdt = (b*(I[i-1]+n*A[i-1])*(S[i-1]/N[i-1])-(s+m+m0) *
                E[i-1])*dt + sd*np.sqrt(dt)*np.random.normal(0, 1)*E[i-1]
        dIdt = (s*(1-w)*E[i-1]-(m+m1+y1)*I[i-1])*dt + \
            sd*np.sqrt(dt)*np.random.normal(0, 1)*I[i-1]
        dAdt = (s*w*E[i-1] - (m+m2+y2)*A[i-1])*dt + \
            sd * np.sqrt(dt)*np.random.normal(0, 1)*A[i-1]
        dRdt = (y1*I[i-1]+y2*A[i-1]-m*R[i-1])*dt + \
            sd*np.sqrt(dt)*np.random.normal(0, 1)*R[i-1]

        S[i] = S[i-1] + dSdt
        E[i] = E[i-1] + dEdt
        I[i] = I[i-1] + dIdt
        A[i] = A[i-1] + dAdt
        R[i] = R[i-1] + dRdt
        N[i] = S[i] + E[i] + I[i] + A[i] + R[i]
    paths.append(E)

# for path in paths:
#     plt.plot(t, path, label='E Something')
# plt.show()

fig, ax = plt.subplots()
ax.plot(t, E, label='Exposed')
ax.plot(t, I, label='Infectious')
ax.plot(t, A, label='Asymptomatic infected')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Number of individuals')
ax.set_title('SEIAR Model')
ax.legend()
plt.show()


# # Plot the results
# fig, ax = plt.subplots()
# ax.plot(t, S, label='Susceptible')
# ax.plot(t, R, label='Recovered')
# ax.set_xlabel('Time (days)')
# ax.set_ylabel('Number of individuals')
# ax.set_title('SEIAR Model')
# ax.legend()
# plt.show()
