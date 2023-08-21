from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

MHz_unit = 1e-3 
GHz_unit = 1

psi01 = tensor(fock(2,1),fock(2,0))
psi02 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 200.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]
w1 = 3.2324 * GHz_unit
w2 = 3.2945 * GHz_unit
wr = 8.2855 * GHz_unit
delta1 = w1 - wr
delta2 = w2 - wr
g1 = 5 * GHz_unit
g2 = g1
g = g1*g2*(delta1+delta2)/(2*delta1*delta2)
alpha = 3.0 * MHz_unit
delta12 = w1 - w2
driving_freq = 2 * MHz_unit
#calculates the hamiltonian of the 1st qubit
H1 = Hamiltonian(1,driving_freq, g, alpha, delta12) 
H2 = Hamiltonian(2,driving_freq, g, alpha, delta12) 

result_1 = mesolve(H1, psi01, times, [], e_ops[1])
result_2 = mesolve(H1,psi02,times,[],e_ops[1])
fig, ax = plt.subplots()
ax.plot(times, result_1.expect[0],label = "|1>")
ax.plot(times,result_2.expect[0],label = "|0>")
ax.set_xlabel('Time')
ax.set_ylabel('Expectation values')
plt.legend()

plt.show()