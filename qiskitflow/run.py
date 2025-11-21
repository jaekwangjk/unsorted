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


#####
##### Prepare physical GRID
#####

# create spatial grid 
L = 30
x = np.linspace(-L/2, L/2, num=nx, endpoint=False)
dx = L / nx

# Gaussian
#sigma = 1
#mu = -4
#g = 1/sigma/np.sqrt(2.*np.pi) * np.exp(-np.square(x-mu)/(2*np.square(sigma)))

sigma = np.sqrt(1/2)
mu = -10 # location of center .... 
g = np.exp(-np.square(x-mu))

for i, v in enumerate(g):
    if np.abs(x[i] - mu) > 5*sigma:
        g[i] = 0

# solve for phi
# initial profile  
phi_init = g
#phi_init_norm = np.sqrt(np.sum(np.square(phi_init)))
#phi0 = phi_init / phi_init_norm


#####
##### Prepare p GRID
#####

# 평행이동, 감쇠 계수, 경계 영역, 주기 폭, why v is used? 
# y 가 p와 같다. 

# v is the warped phase variable 
shift = 0 # shift plays nothing here
alpha = 10 # increase decaying factor for highfrequency component
bzone = 0
wy = 8 # y축방향 길이 
# y>0 뿐만 아니라 y< 0 까지 해서 만든다. (수치적 안정성과 FFT 주기성을 위한 대칭 확장)
# excluse the end point, since sin(y_min)  = sin(y_max) 일것이기때문에 end point를 넣지 않는다. 

# 일단 이해하려는 시도. 
# y가 p에 대응되면
# v는 w가 된다?

# focus on how it carefully treat y (or p) < 0 area 

y = np.linspace(-np.pi*(wy/2+shift), np.pi*(wy/2-shift), num=ny, endpoint=False) # phase 변수의 격자를 만듬. 

# 수치적 안정성을 확보하기 위함. p<0 쪽은 단지 수치적 대칭 확장용이니 지나치게 많은 에너지를 쏟지는 말자. 
bottom = np.exp(-np.pi*(wy/2-shift))

v_init = np.zeros(ny)
for i, vy in enumerate(y):
    if vy < -np.pi*(wy/2+shift-bzone) or vy > np.pi*(wy/2-shift-bzone):
        v_init[i] = 0
    # 비대칭 (convection) -- 단방향적 수송(advective directionality)    
    elif vy < - np.pi * shift :
        tmp = np.exp(alpha*(vy+shift*(1+1/alpha)*np.pi))
        if tmp < bottom :
            tmp = bottom
        v_init[i] = tmp
    else:
        v_init[i] = np.exp(-vy)

v_init0 = np.zeros(ny)
v_init0[:int(ny/2)] = np.exp(y[:int(ny/2)])
v_init0[int(ny/2):] = np.exp(-y[int(ny/2):])


H1, H2= compute_H1_H2_for_Convection_and_Diffusion_and_Source(nx, dx, C ,D, S)

fv = fft.fft(v_init, norm='forward') # 각 k에 해당하는 복소 진폭 및 위상
kv = fft.fftfreq(ny)*ny*2/wy  # fv에 대응하는 주파수 축을 만들어줌, Array of length n containing the sample frequencies.

# prerequisite: H1, H2, kv already computed
M_list = [1j * k * H1 - H2 for k in kv]

# optional: fv 위상 미리 캐싱
phi_phase = np.ones_like(fv, dtype=complex)
nz = np.abs(fv) > 0
phi_phase[nz] = fv[nz] / np.abs(fv[nz])

# t_new=1
# expm(M_list[i]*t_new)

states = []


for i, k in enumerate(kv):
    
    M = 1.j*k*H1 - H2
    U = expm(M*t)
    
    
    
sim_t05=0.5
sim_t10=1.0
sim_t20=2.0


soln_t05 = simulate(phi_init, H1,H2, sim_t05, nx, kv, fv, ny, wy, y,  idx=10)
soln_t10 = simulate(phi_init, H1,H2, sim_t10, nx, kv, fv, ny, wy, y,  idx=10)
soln_t20 = simulate(phi_init, H1,H2, sim_t20, nx, kv, fv, ny, wy, y,  idx=10)


plt.plot(x, soln_t05.real, label="t=0.5")
plt.plot(x, soln_t10.real, label="t=1.0")
plt.plot(x, soln_t20.real, label="t=2.0")

plt.legend()
plt.show()
    