from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import Hamiltonian

psi01 = tensor(fock(2,1),fock(2,0))
psi02 = tensor(fock(2,0),fock(2,0))

print(psi01)
print(psi02)