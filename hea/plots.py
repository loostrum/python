from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np

d='data'
e_bins = np.loadtxt(d+'/bins.txt')
n_sync = np.loadtxt(d+'/sync-input.txt')
n_planck = np.loadtxt(d+'/planck-input.txt')
n_non = np.loadtxt(d+'/non-scattered.txt')
n_comp = np.loadtxt(d+'/scattered.txt')

full_input = np.add(n_planck,n_sync)
full_output = np.add(n_non,n_comp)

plt.loglog(e_bins,full_input,color='blue',label='Sync + planck input')
plt.loglog(e_bins,n_non,color='green',ls='--',label='Non-scattered output')
plt.loglog(e_bins,n_comp,color='red',ls='--',label='Comptonized output')
plt.loglog(e_bins,full_output,color='black',label='Total output')
plt.legend()
plt.xlim(1E-2,1E4)
plt.ylim(1E-35,1E-30)
plt.xlabel(r'E (keV)')
plt.ylabel(r'$\nu \, F_\nu$' )
plt.show()
