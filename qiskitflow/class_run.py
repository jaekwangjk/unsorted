from library import * 


print("PROGRAM EXECUTED")


C = 10
D = 0.5
S = -1

# timestep
t = 1.5

# It solves two-dimensional equation
# number of qubits for u
nqu = 8
q_u = QuantumRegister(nqu, name='qu')
nx = 2**nqu

# number of qubits for v
nqv = 8
ny = 2**nqv


sim = Schrodingeriser(wy=8, ny=2**nqv) # prepare p space


#####
##### Prepare physical GRID


# create spatial grid 
L = 30
x = np.linspace(-L/2, L/2, num=nx, endpoint=False)
dx = L / nx

# Gaussian
sigma = np.sqrt(1/2)
mu = -10 # location of center .... 
g = np.exp(-np.square(x-mu))

for i, v in enumerate(g):
    if np.abs(x[i] - mu) > 5*sigma:
        g[i] = 0

phi_init = g

#### Compute H1 and H2 of finite difference scheme 

H1, H2= compute_H1_H2_for_Convection_and_Diffusion_and_Source(nx, dx, C ,D, S)

soln_t05 = sim.simulate(0.5, phi_init, H1, H2, nx)
soln_t10 = sim.simulate(1.0, phi_init, H1, H2, nx)
soln_t20 = sim.simulate(2.0, phi_init, H1, H2, nx)

plt.plot(x, soln_t05.real, label="t=0.5")
plt.plot(x, soln_t10.real, label="t=1.0")
plt.plot(x, soln_t20.real, label="t=2.0")

plt.legend()
plt.show()
    
    

    