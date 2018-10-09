import numpy as np #this loads a library for linear algebra and comman math functions
import matplotlib.pyplot as plt # this load library for ploting
from numpy import linalg as LA # this loads a library which is in the numpy for linearalgebra this is done for convience
from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#%matplotlib inline
#import mpld3
#mpld3.enable_notebook()
def pm(i): # this is function takes an imput from 0 to 3 and gives an out of pauli matrices such that pm(1)=\sigma_x
    if i==0:
        result=np.matrix([[1., 0.], [0., 1.]]) # this is identity matrix
    elif i==1:
        result=np.matrix([[0., 1.], [1., 0.]]) # sigma_x
    elif i==2:
        result=np.matrix([[0., -1j], [1j, 0.]]) # sigma_y
    elif i==3:
        result=np.matrix([[1., 0.], [0., -1.]]) # sigma_z
    return result
k3=np.arange(-np.pi,np.pi,0.01) # this kreates k3 momentum array from -pi to pi
k2=0. # set k2 zero
k1=0. # and set k1 zero
N=50 # this is number of sites
gamma=0. # m1 and m2 controls the distance between the weyl nodes
lam=-1.
tx=1.
t=1.
k0=np.pi/3
#for i, val in enumerate(k1):
#    for j, val in enumerate(k3):
#        H=-(m*(2.-np.cos(k2)-np.cos(k3[j]))+2*tx*(np.cos(k1[i])-np.cos(k0)))*pm(1)-2*t*np.sin(k2)*pm(2)-2*t*np.sin(k3[j])*pm(3)
#        E=LA.eigvalsh(H)
#        ax.plot(np.ones((np.size(E)))*k1[i],np.ones((np.size(E)))*k3[j],E,'ko',markersize=0.5)
fig1=plt.figure(1)
for i, val in enumerate(k3):
    H=(lam*(2.-np.cos(k1)-np.cos(k2))+(np.cos(k3[i])-np.cos(k0)))*pm(3)+np.sin(k2)*pm(2)+np.sin(k1)*pm(1)
    E=LA.eigvalsh(H)
    plt.plot(np.ones((np.size(E)))*k3[i],E,'ko',markersize=0.5)

fig2=plt.figure(2)
#super conducting case
for i, val in enumerate(k3):
    H=(lam*(2.-np.cos(k1)-np.cos(k2))+(np.cos(k3[i])-np.cos(k0)))*pm(3)+np.sin(k2)*pm(2)+np.sin(k1)*pm(1)
    Hm=(lam*(2.-np.cos(-k1)-np.cos(-k2))+(np.cos(-k3[i])-np.cos(k0)))*pm(3)+np.sin(-k2)*pm(2)+np.sin(-k1)*pm(1)
    h_up=np.hstack((H,np.zeros((2,2))))    # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
    h_down=np.hstack((np.zeros((2,2)),-np.transpose(Hm)))
    h_bog=np.vstack((h_up,h_down))
    delta=0.3
    slope=1
    delta_up=np.hstack((np.zeros((2,2)),delta*1j*pm(2)/slope))    # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
    delta_down=np.hstack((np.conj(delta)*-1j*pm(2)/slope,np.zeros((2,2))))
    delta=np.vstack((delta_up,delta_down))
    h_bog=delta+h_bog
    E=LA.eigvalsh(h_bog)
    plt.plot(np.ones((np.size(E)))*k3[i],E,'ko',markersize=0.5)
plt.show()
print(delta)
