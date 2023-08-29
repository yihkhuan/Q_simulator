from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian,rabi, Hamiltonian_phase


MHz_unit = 1e-3 
GHz_unit = 1

psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 1000.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]
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
eta = 0 
phi = 0
driving_x = 2 *MHz_unit
#calculates the hamiltonian of the 1st qubit
H1 = Hamiltonian(1,driving_freq, g, alpha, delta12, eta, driving_x, phi) 
H2 = Hamiltonian(2,driving_freq, g, alpha, delta12, eta, driving_x, phi) 

result_gnd = mesolve(H1, psi01, times, [], e_ops)
#result_2 = mesolve(H1, psi02, times, [], e_ops[0])
#result_3 = mesolve(H1, psi01, times, [], e_ops[1])
result_exc = mesolve(H1, psi02, times, [], e_ops)

rabi_gnd = rabi(times,result_gnd.expect[1])
rabi_exc = rabi(times,result_exc.expect[1])

print(rabi_gnd, end = ", ")
print(rabi_exc, end = ", ")
print(rabi_gnd - rabi_exc)



# plt.figure()
# plt.plot(times, result_gnd.expect[0],label = "Q1 = |0>")
# plt.plot(times,result_exc.expect[0],label = "Q1 = |1>")
# plt.xlabel('Time, ns')
# plt.ylabel('Expectation values')
# plt.title("Q1")
# plt.legend()

plt.figure()
plt.plot(times, result_gnd.expect[1],label = "Q1 = |0>")
plt.plot(times,result_exc.expect[1],label = "Q1 = |1>")
plt.xlabel('Time, ns')
plt.ylabel('Expectation values')
plt.title("wr = 7.5GHz, g1 = 20 MHz")
plt.legend()


plt.show()

