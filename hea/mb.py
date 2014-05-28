import scipy.stats as stats
from matplotlib import pyplot as plt
import numpy as np

maxwell=stats.maxwell
# loc = shift from 0, scale = a on wikipedia, size is amount of numbers
data=maxwell.rvs(loc=0, scale=5, size=1E6)
plt.hist(data,bins=100)
x=np.linspace(0,30,100)
plt.show()
