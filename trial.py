from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

# psi01 = tensor(fock(2,1),fock(2,0))
# psi02 = tensor(fock(2,0),fock(2,0))

# times = np.linspace(0.0,10.0,1000)

# psi = (2.0 * basis(2, 0) + basis(2, 1)).unit()

# def linear(t, args):
#     return 0.3 * t

# def periodic(t, args):
#     return np.cos(0.5 * t)

# # Define QobjEvos
# H_lin = QobjEvo([[sigmaz(), linear]], tlist=times)
# H_per = QobjEvo([[sigmaz(), periodic]], tlist=times)

# H = [H_lin,H_per]

# result_lin = sesolve(H, psi, times, [sigmay()])

# # Plot <sigma_y> for linear increasing field strength
# plt.plot(times, result_lin.expect[0])
# plt.xlabel("Time"), plt.ylabel("<sigma_y>")

# plt.plot(times, result_lin.expect[1])
# plt.xlabel("Time"), plt.ylabel("<sigma_y>")
# plt.show()
xlist = np.linspace(0, 4.0, 5)
ylist = np.linspace(8.0, 12.0, 5)
X, Y = np.meshgrid(xlist, ylist)
Z = (X+Y)
print(Z)


