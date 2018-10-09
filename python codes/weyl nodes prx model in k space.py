import numpy as np  # this loads a library for linear algebra and comman math functions
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import \
    linalg as LA  # this loads a library which is in the numpy for linearalgebra this is done for convience


# %matplotlib inline
# import mpld3
# mpld3.enable_notebook()
def pm(i):  # this is function takes an imput from 0 to 3 and gives an out of pauli matrices such that pm(1)=\sigma_x
    if i == 0:
        result = np.matrix([[1., 0.], [0., 1.]])  # this is identity matrix
    elif i == 1:
        result = np.matrix([[0., 1.], [1., 0.]])  # sigma_x
    elif i == 2:
        result = np.matrix([[0., -1j], [1j, 0.]])  # sigma_y
    elif i == 3:
        result = np.matrix([[1., 0.], [0., -1.]])  # sigma_z
    return result


def tridiag(h1, h2, h3,
            N):  # this creates a tridiagonal block matrix main diagonal consists of h1 first diag above is h2 and first diag below is h3
    h = np.kron(np.identity(N), h1)  # this creates a block diagonal matrix with h1 is diagonal elements
    spm = np.diag(np.ones(N - 1), 1)
    hp = np.kron(spm, h2)  # this replaces first above diag elements with h2's
    spm = np.diag(np.ones(N - 1), -1)
    hm = np.kron(spm, h3)  # this replaces first below diag elements with h3's
    result = h + hm + hp  # this creates the desired matrix
    return result


def kronsum(a, b):
    a_up = np.hstack((a, np.zeros((np.size(a, 0), np.size(b,
                                                          1)))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
    b_down = np.hstack((np.zeros((np.size(b, 0), np.size(a, 1))), b))
    result = np.vstack((a_up, b_down))
    return result


def Mf(i, n, N,range):
    if n == 1:
        result = range #2.15 is good
    elif n == 2:
        result = i
    elif n == 3:
        result = np.exp(-((i - N) / 50) ** 2)
    elif n == 4:
        if (i < N / 2):
            result = range*i/(N/2)
        else:
            result = range*(N - i - 1)/(N/2)

    return result


k3 = np.arange(-np.pi, np.pi, 0.1)  # this kreates k1 momentum array from -pi to pi
k1 = k3  # set k2 zero
k2 = 0.  # and set k3 zero
N = 50  # this is number of sites

M= 1.15

t = 1.
fs=1
rr=2.15

fig = plt.figure()
ax = fig.gca(projection='3d')
for i, val in enumerate(k3):
    for j, val in enumerate(k1):
        h = -pm(2)*np.sin(k2)-pm(3)*(np.cos(k1[j])+np.cos(k2)+np.cos(k3[i])-M)-pm(1)*np.sin(k1[j])

        w = LA.eigvalsh(h)  # this calculates the eigen values aka energies of the system for given k_3

        ax.plot(np.ones((np.size(w))) * k3[i]/np.pi, np.ones((np.size(w))) * k1[j]/np.pi,w, 'ko', markersize=0.1)  # this puts dots for each eigen value for a given k_3





plt.show()