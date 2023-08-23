from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

def mu(sign: int, g:float, alpha: float, delta12: float):
    #mu = sign*(g*alpha/(delta12*(alpha-sign*delta12)))
    mu = sign * g / delta12 * alpha / (alpha - sign * delta12)
    return mu

def nu(sign: int, g:float, alpha: float, delta12: float):
    #nu = -(g/(alpha-sign*delta12))
    nu = sign * g / delta12 * (-sign * delta12) / (alpha - sign * delta12)
    return nu

def Hamiltonian (QB: int, driving: float, g: float, alpha: float, delta12: float,eta):
    H = np.zeros(tensor(identity(2),sigmax()).shape[0])
    match QB:
        case 1:
            H = driving*(nu(-1,g,alpha,delta12)*tensor(identity(2),sigmax())
                         +mu(-1,g,alpha,delta12)*tensor(sigmaz(),sigmax())
                         +eta*tensor(sigmaz(),identity(2)))
        case 2:
            H = driving*(nu(1,g,alpha,delta12)*tensor(sigmax(),identity(2))
                         +mu(1,g,alpha,delta12)*tensor(sigmax(),sigmaz())
                         +eta*tensor(identity(2),sigmaz()))
    return H

def D_Hamiltonian(driving: float, phi: float):
    H = -driving*(np.cos(phi)*sigmax()+np.sin(phi)*sigmay())
    return H

#function outputs the rabi frequency of the qubit
def rabi (times: np.ndarray, plot: list) -> float:
    N = np.size(times)
    T = times[1]
    yf = fft(plot)
    xf = fftfreq(N, T)[:N//2]
    nmax = np.argmax(np.abs(yf[0:N//2]))
    rb_freq = xf[nmax]
    return rb_freq