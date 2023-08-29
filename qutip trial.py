from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import D_Hamiltonian

psi0 = fock(2,1)
#time list to measure the qubit
times = np.linspace(0.0, 10.0, 1000)
#expectation operator for the 1st qubit
e_ops = sigmaz()
driving_freq = 5.0 
g = 2.0
alpha = 3.0
delta12 = 1.0
phi = np.pi/3.0
#calculates the hamiltonian of the 1st qubit
H = D_Hamiltonian(5.0,phi)

result = mesolve(H, psi0, times, [], e_ops)
plt.subplots()
plt.plot(times, result.expect[0])
plt.set_xlabel('Time')
plt.set_ylabel('Expectation values')
plt.legend(("Sigma-Z"))
plt.show()
