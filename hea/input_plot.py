import numpy as np
from matplotlib import pyplot as plt

data_sync = np.loadtxt('sync.txt')
sync_bins = data_sync[0]
sync_weights = data_sync[1]

data_planck = np.loadtxt('planck.txt')
planck_bins = data_planck[0]
planck_weights = data_planck[1]

plt.loglog(sync_bins,sync_weights,label='sync')
plt.loglog(planck_bins,planck_weights,label='planck')
plt.legend(loc=4)
plt.ylim(1E-10,1E30)
plt.show()
