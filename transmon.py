from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian_transmon

MHz_unit = 1e-3 
GHz_unit = 1

psi02 = fock(2,0)
psi01 = fock(2,1)
superpsi = (1*fock(2,0)+np.sqrt(2)*fock(2,1)) / np.sqrt(3)
#time list to measure the qubit
times = np.linspace(0.0, 1000.0, 10000)
#expectation operator for the 1st qubit
e_ops = sigmaz()
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
driving_freq =  20000 * MHz_unit

#Hamiltonian transmon
H = Hamiltonian_transmon(w1,w2,alpha,g)

#solve SE
results01 = mesolve(H,psi01,times,[],e_ops)
results02 = mesolve(H,psi02,times,[],e_ops)
resultssuper = mesolve(H,superpsi,times,[],e_ops)

#plot results
plt.figure()
plt.plot(times,results01.expect[0],label = "Qubit 1")
plt.plot(times,results02.expect[0],label = "Qubit 2")
plt.plot(times,resultssuper.expect[0],label = "Qubit 3")
plt.xlabel('time, ns')
plt.ylabel('Expectation value')
plt.title('Transmon qubit hamiltonian')
plt.legend()
plt.grid()
plt.show()

