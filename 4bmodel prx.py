import numpy as np  # this loads a library for linear algebra and comman math functions
import matplotlib.pyplot as plt  # this load library for ploting
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
def g(i):
    if i==1:
        result=np.kron(pm(3),pm(2))
    elif i==2:
        result = np.kron(pm(3), pm(1))
    elif i == 3:
        result = np.kron(pm(2), pm(0))
    elif i == 4:
        result = np.kron(pm(1), pm(0))
    elif i == 5:
        result = np.kron(pm(2), pm(3))
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
            result = range*(N - i - 1.)/(N/2)

    return result


k1 = np.arange(-0.4*np.pi, 0.4*np.pi, 0.01)  # this kreates k1 momentum array from -pi to pi
k2 = 0.0  # set k2 zero
k3 = 0.0  # and set k3 zero
N = 100  # this is number of sites
m=.5
M= 3.+m

t = 1.
fs=4
rr=1.

fig1 = plt.figure(1)
plt.subplot(251)
for i, val in enumerate(k1):
    h1 = -g(1)*np.sin(k1[i])-g(3)*np.sin(k3)-(np.cos(k1[i])+np.cos(k3)-M)*g(4)
    h2 = -g(4) / 2. + g(2) / 2.*1j
    h3 = -g(4) / 2. - g(2) / 2.*1j
    H= tridiag(h1, h2, h3, N)  # this creates 2N*2N hamiltonian that describes the system for a GIVEN k_3
    b = np.empty(shape=(0, 0))
    for ii in range(N):
        b=kronsum(b,g(5)*g(1)*Mf(ii,fs,N,rr))
    H=H+b
    w = LA.eigvalsh(H)  # this calculates the eigen values aka energies of the system for given k_3

    plt.plot(np.ones((np.size(w))) * k1[i]/np.pi, w, 'ko',
             markersize=0.1)  # this puts dots for each eigen value for a given k_3
    plt.xlim((-0.6, 0.6))
    plt.ylim((-1.2, 1.2))

fig2=plt.figure(2)
plt.subplot(252)
for ii in range(N):
    plt.plot(ii,Mf(ii,fs,N,rr),'ko',markersize=3)


plt.show()