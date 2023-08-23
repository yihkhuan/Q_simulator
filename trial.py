from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))

times = np.linspace(0.0,100.0,10000)

# Define QobjEvos
H = tensor(qeye(2),sigmax())
e_ops = tensor(qeye(2),sigmaz())

result_gnd = mesolve(H, psi01, times, [], e_ops)
result_exc = mesolve(H, psi02, times, [], e_ops)

plt.figure()
plt.plot(times, result_gnd.expect[0],label = "Q1 = |0>")
plt.plot(times,result_exc.expect[0],label = "Q1 = |1>")
plt.xlabel('Time, ns')
plt.ylabel('Expectation values')
plt.title("Q1")
plt.legend()
plt.show()


