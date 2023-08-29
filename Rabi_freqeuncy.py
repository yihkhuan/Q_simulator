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

N = 15
w1 = 4.3 * GHz_unit
w2 = 4.5 * GHz_unit
wr = np.linspace(6.0,10.0,N) * GHz_unit
g1 = np.linspace(10.0,100.0,N) * MHz_unit
alpha = -220 * MHz_unit
delta12 = w1 - w2
driving_freq =  20000 * MHz_unit
eta = 0
#calculates the hamiltonian of the ith qubit and the resulting output of the hamiltonian
X, Y = np.meshgrid(wr, g1)
results = np.zeros(X.shape)

i = 0

for glist in g1:
    j = 0
    for w in wr:
        delta1 = w1 - w
        delta2 = w2 - w
        g = (glist*glist)*(delta1+delta2)/(2*delta1*delta2)
        H = Hamiltonian(1,driving_freq, g, alpha, delta12, eta)
        result_gnd = mesolve(H, psi01, times, [], e_ops)
        result_exc = mesolve(H, psi02, times, [], e_ops)
        rabi_gnd = rabi(times,result_gnd.expect[1])
        rabi_exc = rabi(times,result_exc.expect[1])
        deltarabi = rabi_gnd - rabi_exc
        results[i][j] = deltarabi
        j = j + 1
    i = i + 1

nslice = N + 1
a = np.arange(nslice)
b = np.max(results) / (nslice-1)

#plotting
plt.subplots(1,1)
cp = plt.contourf(X, Y, results,b*a)
plt.colorbar(cp,label = 'Rabi difference, GHz') # Add a colorbar to a plot
plt.title('Delta Rabi plot')
plt.xlabel('resonance frequency, GHz')
plt.ylabel('g1, GHz')
plt.show()

# delta1 = w1 - X
# delta2 = w2 - X
# g = (Y*Y)*(delta1+delta2)/(2*delta1*delta2)

# Rabi_gnd = driving_freq*g/delta12
# Rabi_exc = driving_freq*g/(delta12*(alpha+delta12))*(delta12-alpha)

# deltarabi = Rabi_gnd - Rabi_exc

# nslice = N + 1
# a = np.arange(nslice)
# b = np.max(deltarabi) / (nslice-1)

# plt.subplots(1,1)
# cp = plt.contourf(X, Y, deltarabi, a*b)
# plt.colorbar(cp,label = 'Rabi difference, GHz') # Add a colorbar to a plot
# plt.title('Delta Rabi plot')
# plt.xlabel('resonance frequency, GHz')
# plt.ylabel('g1, GHz')
# plt.show()

