from scipy.fft import fft, fftfreq
import numpy as np

#function outputs the rabi frequency of the qubit
def rabi (times: np.ndarray, plot: list) -> float:
    N = np.size(times)
    T = 1.0/times[1]
    yf = fft(plot)
    xf = fftfreq(N, T)[:N//2]
    nmax = np.argmax(np.abs(yf[0:N//2]))
    rb_freq = xf[nmax]
    return rb_freq

#Qubit simulation
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

MHz_unit = 1e-3 
GHz_unit = 1
#initial state of qubits and the parameters
psi02 = tensor(fock(2,1),fock(2,0))
psi01 = tensor(fock(2,0),fock(2,0))
#time list to measure the qubit
times = np.linspace(0.0, 1800.0, 10000)
#expectation operator for the 1st qubit
e_ops = [tensor(sigmaz(),identity(2)),tensor(identity(2),sigmaz())]

N = 50
w1 = 4.3 * GHz_unit
w2 = 4.5 * GHz_unit
wr = np.linspace(6.0,10.0,N) * GHz_unit
delta1 = w1 - wr
delta2 = w2 - wr
g1 = np.linspace(10.0,100.0,N) * MHz_unit
g2 = g1
glist = g1*g2*(delta1+delta2)/(2*delta1*delta2)
alpha = -220 * MHz_unit
delta12 = w1 - w2
driving_freq =  20000 * MHz_unit
eta = 0
#calculates the hamiltonian of the ith qubit and the resulting output of the hamiltonian
X, Y = np.meshgrid(wr, g1)
results = np.zeros(X.shape)
i = 0

for w in wr:
    j = 0
    for g in glist:
        H = Hamiltonian(1,driving_freq, g, alpha, delta12, eta)
        result_gnd = mesolve(H, psi01, times, [], e_ops)
        result_exc = mesolve(H, psi02, times, [], e_ops)
        rabi_gnd = rabi(times,result_gnd.expect[1])
        rabi_exc = rabi(times,result_exc.expect[1])
        deltarabi = rabi_gnd - rabi_exc
        results[i][j] = deltarabi
        print(i,end = ",")
        print(j)

        j = j + 1
    i = i + 1

#plotting
plt.subplots(1,1)
cp = plt.contourf(X, Y, results)
plt.colorbar(cp,label = 'Rabi difference, GHz') # Add a colorbar to a plot
plt.title('Delta Rabi plot')
plt.xlabel('resonance frequenxy, GHz')
plt.ylabel('g1, MHz')
plt.show()


