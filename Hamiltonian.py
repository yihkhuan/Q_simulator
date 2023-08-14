from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def mu(sign: int, g:float, alpha: float, delta12: float):
    mu = sign*(g*alpha/(delta12*(alpha-sign*delta12)))
    return mu

def nu(sign: int, g:float, alpha: float, delta12: float):
    nu = sign*(g*(-sign)*delta12/(delta12*(alpha-sign*delta12)))
    return nu

def Hamiltonian (QB: int, driving: float, g: float, alpha: float, delta12: float):
    H = np.zeros(tensor(identity(2),sigmax()).shape[0])
    match QB:
        case 1:
            H = driving*(tensor(sigmax(),identity(2))+nu(-1,g,alpha,delta12)*tensor(identity(2),sigmax())
                         +mu(-1,g,alpha,delta12)*tensor(sigmaz(),sigmax()))
        case 2:
            H = driving*(tensor(identity(2),sigmax())+nu(1,g,alpha,delta12)*tensor(sigmax(),identity(2))
                         +mu(1,g,alpha,delta12)*tensor(sigmax(),sigmaz()))
    return H

def D_Hamiltonian(driving: float, phi: float):
    H = -driving*(np.cos(phi)*sigmax()+np.sin(phi)*sigmay())
    return H