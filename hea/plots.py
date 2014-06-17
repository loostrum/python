from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np


e_bins = np.loadtxt('bins.txt')
n_sync = np.loadtxt('sync-input.txt')
n_planck = np.loadtxt('planck-input.txt')
n_non = np.loadtxt('non-scattered.txt')
n_comp = np.loadtxt('scattered.txt')

full_input = np.add(n_planck,n_sync)
full_output = np.add(n_non,n_comp)

plt.loglog(e_bins,full_input,color='blue',ls='--',label='Sync + planck input')
plt.loglog(e_bins,full_output,color='red',label='Comptonized output')
plt.legend()
plt.xlim(1E-2,1E4)
plt.ylim(1E-35,1E-30)
plt.xlabel(r'E (keV)')
plt.ylabel(r'$\nu \, F_\nu$' )
plt.show()
