import numpy as np

from scipy import fft, special
from scipy.linalg import expm, eig
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit_aer import Aer  # 또는 AerSimulator
from qiskit.compiler import transpile
from qiskit.quantum_info import Operator, state_fidelity
from qiskit.circuit.library import QFT

class PlotConfig:
    def __init__(self, figsize=(6,4), dpi=120, subplot_ratio=None):
        self.figsize = figsize
        self.dpi = dpi
        self.subplot_ratio = subplot_ratio

    def get_simple(self):
        if self.subplot_ratio:  
            # subplot_ratio 활용해서 가로/세로 비율 조정
            fig, ax = plt.subplots(
                figsize=(self.figsize[0], self.figsize[1]*self.subplot_ratio),
                dpi=self.dpi
            )
        else:
            fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        ax.grid(True, alpha=0.3)
        return fig, ax
   

def analytic_solution(t):
    return np.exp(-(x - mu - C*t)**2 / (1 + 4*D*t)) \
           / np.sqrt(1 + 4*D*t) \
           * np.exp(S*t)
           


def compute_H1_H2_for_Convection_and_Diffusion_and_Source(nx, dx, C ,D, S):

    AC = np.zeros((nx, nx))
    
    for i in range(nx-1):
        AC[i,i+1] = 1
        AC[i+1,i] = -1
    
    AC[0,-1] = -1
    AC[-1,0] = 1    
    
    AC = AC * C / (2*dx)
    
    AD = np.zeros((nx, nx))
    for i in range(nx-1):
        AD[i,i+1] = -1
        AD[i+1,i] = -1
        AD[i,i] = 2
    AD[-1,-1] = 2
    AD[0,-1] = -1
    AD[-1,0] = -1
    
    AD = AD * D / np.square(dx)
    
    # Source Term

    AS = np.eye(nx)
    AS *= -S


    # Construct overall operator 
    A = AC + AD + AS

    A = np.matrix(A)
    AH = A.H # Hermitian conjugate (즉, 전치 + 켤레)
    H1 = (A+AH)/2  # Hermitian (자기수반) 부분
    H2 = (A-AH)/2 # Skew-Hermitian (비자기수반) 부분
    
    return H1, H2
    
    
def simulate(phi_init, H1,H2, t, nx, kv, fv, ny,wy,y, idx=10):
    
    
    phi_init_norm = np.sqrt(np.sum(np.square(phi_init)))
    phi0 = phi_init / phi_init_norm
    
    
    for i, k in enumerate(kv):
    
        M = 1.j*k*H1 - H2
        U = expm(M*t)
    
    # y-슬라이스 인덱스
    y_idx = int(ny/2 + idx)

    soln = np.zeros(nx, dtype='complex128')
    

    for i, k in enumerate(kv):

        # ---- (1) k-모드 초기화 (양자 initialize 대응) ----
        if np.abs(fv[i]) > 0:
            phase = fv[i] / np.abs(fv[i])
        else:
            phase = 1.0 + 0j
        init_k = phi0 * phase                      # shape = (nx,)

        # ---- (2) 시간 진화 U(t) ----
        M = 1j * k * H1 - H2                       # (nx x nx)
        U = expm(M * t)                            # (nx x nx)
        state_t = U @ init_k                       # evolved mode

        # ---- (3) k-모드 합성 (inverse transform 기준) ----
        soln += state_t * np.abs(fv[i]) * np.exp(
            1j * np.pi * k * wy/2 * 2 * y_idx / ny
        )

    # ---- (4) 초기 정규화 복원 ----
    soln *= phi_init_norm
    soln *= np.exp(y[y_idx])      # (원래 코드 그대로 반영)

    return soln