import numpy as np
import matplotlib.pyplot as plt  # this load library for ploting
ef=np.arange(0.1,0.3,0.001)
wd=0.2
g0=350.
k0=np.pi/3
vf=np.sin(k0)
delta=5*wd*np.exp(-6*np.pi**2*vf*2/ef**2/g0)
plt.plot(ef,delta)
plt.show()