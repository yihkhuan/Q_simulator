from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian_RWA,Hamiltonian_5


MHz_unit = 1e-3 
GHz_unit = 1

psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 1000.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]
w1 = 4.3 * GHz_unit
w2 = 4.65 * GHz_unit
wr = 6.5 * GHz_unit
delta1 = w1 - wr
delta2 = w2 - wr
g1 = 50 * MHz_unit
g2 = g1
g = g1*g2*(delta1+delta2)/(2*delta1*delta2)
g = 10 * MHz_unit
delta12 = w1 - w2
driving_freq =  500 * MHz_unit
alpha = -220 * MHz_unit

#calculates the hamiltonian of the 1st qubit
H1 = Hamiltonian_RWA(driving_freq, g, w1-w2, w2 - w2)
H2 = Hamiltonian_5(w1,w2,delta12,alpha,g,driving_freq)
result_gnd = mesolve(H2, psi01, times, [], e_ops)
#result_2 = mesolve(H1, psi02, times, [], e_ops[0])
#result_3 = mesolve(H1, psi01, times, [], e_ops[1])
result_exc = mesolve(H2, psi02, times, [], e_ops)



plt.figure()
plt.plot(times, result_gnd.expect[0],label = "Q1 = |0>")
plt.plot(times,result_exc.expect[0],label = "Q1 = |1>")
plt.xlabel('Time, ns')
plt.ylabel('Expectation values')
plt.title("Q1")
plt.legend()

plt.figure()
plt.plot(times, result_gnd.expect[1],label = "Q1 = |0>")
plt.plot(times,result_exc.expect[1],label = "Q1 = |1>")
plt.xlabel('Time, ns')
plt.ylabel('Expectation values')
plt.title("Q2")
plt.legend()


plt.show()

