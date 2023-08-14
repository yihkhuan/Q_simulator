from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

MHz_unit = 1E6 
ns_unit = 1E-9

psi01 = tensor(fock(2,1),fock(2,0))
psi02 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 5.0, 1000)
#expectation operator for the 1st qubit
e_ops1 = tensor(sigmaz(),identity(2))
e_ops2 = tensor(identity(2),sigmaz())
driving_freq = 3
alpha = 3.0
delta12 = 0.5
g = 0.0001 * delta12
#calculates the hamiltonian of the 1st qubit
H1 = Hamiltonian(1,driving_freq, g, alpha, delta12) 
H2 = Hamiltonian(2,driving_freq, g, alpha, delta12) 

result_1 = mesolve(H1, psi01, times, [], e_ops1)
result_2 = mesolve(H1,psi02,times,[],e_ops1)
fig, ax = plt.subplots()
ax.plot(times, result_1.expect[0],label = "|1>")
ax.plot(times,result_2.expect[0],label = "|0>")
ax.set_xlabel('Time')
ax.set_ylabel('Expectation values')
plt.legend()

plt.show()