#Qubit simulation
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian,rabi

MHz_unit = 1e-3 
GHz_unit = 1
#initial state of qubits and the parameters
psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 10000.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]

N = 100
w1 = 4.3 * GHz_unit
w2 = 4.5 * GHz_unit
wr = 7.5 * GHz_unit
g1 = 60 * MHz_unit
alpha = np.linspace(0.0,500.0,N) * MHz_unit
delta12 = w1 - w2
driving_freq =  20000 * MHz_unit
eta = 0
#calculates the hamiltonian of the ith qubit and the resulting output of the hamiltonian

results_gnd = np.zeros(N)
results_exc = np.zeros(N)
i = 0

delta1 = w1 - wr
delta2 = w2 - wr
g = (g1*g1)*(delta1+delta2)/(2*delta1*delta2)

for a in alpha:
    H = Hamiltonian(1,driving_freq, g, a, delta12, eta)
    result_gnd = mesolve(H, psi01, times, [], e_ops[1])
    result_exc = mesolve(H, psi02, times, [], e_ops[1])
    rabi_gnd = rabi(times,result_gnd.expect[0])
    rabi_exc = rabi(times,result_exc.expect[0])
    results_gnd [i] = rabi_gnd
    results_exc[i] = rabi_exc
    print(i)
    i = i + 1
    

plt.plot(alpha, results_gnd,label = "ground")
plt.plot(alpha, results_exc,label = "excited")
plt.title('alpha vs rabi(exc)')
plt.xlabel('alpha, GHz')
plt.ylabel('rabi frequency of qubit 2, GHz')
plt.legend()
plt.show()


