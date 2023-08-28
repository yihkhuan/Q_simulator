from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian,rabi
from scipy.fft import fft, fftfreq


MHz_unit = 1e-3 
GHz_unit = 1

psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 10000.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]
w1 = 4.3 * GHz_unit
w2 = 4.5 * GHz_unit
wr = 9 * GHz_unit
delta1 = w1 - wr
delta2 = w2 - wr
g1 = 80 * MHz_unit
g2 = g1
g = g1*g2*(delta1+delta2)/(2*delta1*delta2)
alpha = -220 * MHz_unit
delta12 = w1 - w2
driving_freq =  20000 * MHz_unit
eta = 0.2 
#calculates the hamiltonian of the 1st qubit
H1 = Hamiltonian(1,driving_freq, g, alpha, delta12, eta,0,0) 
H2 = Hamiltonian(2,driving_freq, g, alpha, delta12, eta,0,0) 

result_gnd = mesolve(H1, psi01, times, [], e_ops)
#result_2 = mesolve(H1, psi02, times, [], e_ops[0])
#result_3 = mesolve(H1, psi01, times, [], e_ops[1])
result_exc = mesolve(H1, psi02, times, [], e_ops)

rabi_gnd = rabi(times,result_gnd.expect[1])
rabi_exc = rabi(times,result_exc.expect[1])


plot_gnd = result_gnd.expect[1]
plot_exc = result_exc.expect[1]

N = np.size(times)
T = times[1]
yf_gnd = fft(plot_gnd)
yf_exc = fft(plot_exc)
xf = fftfreq(N, T)[:N//2]
a = np.abs(yf_gnd[0:N//2])
b = np.abs(yf_exc[0:N//2])

plt.plot(xf, 2.0/N * np.abs(yf_gnd[0:N//2]), label = "gnd")
plt.plot(xf, 2.0/N * np.abs(yf_exc[0:N//2]), label = "exc")
plt.legend()
plt.grid()
plt.show()


