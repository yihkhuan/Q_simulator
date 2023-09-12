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

def Hamiltonian (QB: int, driving: float, g: float, alpha: float, delta12: float, eta) -> Qobj:
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

def Hamiltonian_phase(driving : float, phi : float) -> Qobj :
    H = -driving*1/2*(np.cos(phi)*tensor(identity(2),sigmax())+np.sin(phi)*tensor(identity(2),sigmay()))
    return H

def Hamiltonian_transmon(w1:float, w2: float, alpha: float, g: float) -> Qobj :
    b = tensor(fock(2,0),dag(fock(2,1)))
    H = g*(dag(b)*b+b*dag(b))
    ws = [w1,w2]
    for w in ws:
        H = H + w*dag(b)*b+alpha/2*dag(b)*b*(dag(b)*b-1)
    return H

def Hamiltonian_CR(driving:float, delta12:float, g:float) -> Qobj :
    H = ((delta12 - np.sqrt(delta12 ** 2 + driving ** 2)) / 2 * tensor(sigmaz(),qeye(2)) - 
         g*driving / np.sqrt(delta12 ** 2 + driving ** 2) / 2 * tensor(sigmaz(),sigmax()))
    return H
    
def Hamiltonian_RWA(driving:float, g: float, delta1: float, delta2: float) -> Qobj :
    H1 = (delta1) / 2 * tensor(sigmaz(),identity(2)) + g / 2 * (tensor(sigmax(),sigmax()) + tensor(sigmay(),sigmay())) + driving * tensor(identity(2), sigmax()) / 2
    H2 = (delta2) / 2 * tensor(identity(2),sigmaz()) + g / 2 * (tensor(sigmax(),sigmax()) + tensor(sigmay(),sigmay())) + driving * tensor(identity(2), sigmax()) / 2
    H = H1 + H2
    return H

def Hamiltonian_5(w1 ,w2, delta12, alpha, J, driving) -> Qobj:
    wix = -J*driving/(delta12 + alpha)
    wzi = (1/(2 * (delta12 + alpha)) - 1/(2 * delta12))*driving**2
    wzx = (1/(delta12 + alpha) - 1/(delta12)) * J * driving
    wzz = (1/(delta12 - alpha) - 1/(delta12 + alpha)) * J **2

    H = (wix * tensor(identity(2),sigmax()) +
         wzi * tensor(sigmaz(),identity(2)) +
         wzx * tensor(sigmaz(),sigmax()) +
         wzz * tensor(sigmaz(),sigmaz())) / 2
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