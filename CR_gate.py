from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian_CR

MHz_unit = 1e-3 
GHz_unit = 1

psi00 = tensor(fock(2,0),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,1))
psi10 = tensor(fock(2,1),fock(2,0))
psi11 = tensor(fock(2,1),fock(2,1))
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
driving = 500 * GHz_unit

#hamiltonian
H = Hamiltonian_CR(driving,delta12,g)

#result
results00 = mesolve(H,psi00,times,[],e_ops)
results01 = mesolve(H,psi01,times,[],e_ops)
results10 = mesolve(H,psi10,times,[],e_ops)
results11 = mesolve(H,psi11, times, [], e_ops) 
plt.figure
plt.plot(times,results00.expect[0],label = "00")
plt.plot(times,results01.expect[0],label = "01")
plt.plot(times,results10.expect[0],label = "10")
plt.plot(times,results11.expect[0],label = "11")
plt.xlabel('Time,ns')
plt.ylabel('Expectation value')
plt.title('Expectation vs time')
plt.grid()
plt.legend()
plt.show()
