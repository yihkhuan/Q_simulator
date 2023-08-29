from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian_phase

MHz_unit = 1e-3 
GHz_unit = 1

psi0 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 5000.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(identity(2),sigmax()), tensor(identity(2),sigmay()), tensor(identity(2),sigmaz())]
w1 = 4.3 * GHz_unit
w2 = 4.5 * GHz_unit
wr = 6.5 * GHz_unit
delta1 = w1 - wr
delta2 = w2 - wr
g1 = 40 * MHz_unit
g2 = g1
g = g1*g2*(delta1+delta2)/(2*delta1*delta2)
alpha = -220 * MHz_unit
delta12 = w1 - w2
eta = 0 
phi = 0
driving_x = 2 *MHz_unit

#Hamiltonian
H = Hamiltonian_phase(driving_x,phi)

#Solve SE
results = mesolve(H,psi0,times,[],e_ops)

#display results
plt.figure
plt.plot(times,results.expect[0],label = "x")
plt.plot(times,results.expect[1],label = "y")
plt.plot(times,results.expect[2],label = "z")
plt.xlabel('Time,ns')
plt.ylabel('Expectation value')
plt.title('Expectation vs time')
plt.grid()
plt.legend()
plt.show()
