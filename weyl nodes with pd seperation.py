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


def b(i, n, N,range):
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


k3 = np.arange(-np.pi/2., np.pi/2., 0.01/2.)  # this kreates k1 momentum array from -pi to pi
k3=np.array([0.])
k2 = 0.  # set k2 zero
k1 = 0.  # and set k3 zero
N = 70  # this is number of sites


t = 1.
fs=4
rr=5.
mg=1.1
M=3.*0.
fig1 = plt.figure(1)
for i, val in enumerate(k3):
    h1 = -pm(1)*np.sin(k1)-pm(3)*(np.cos(k1)+np.cos(k3[i])-mg-M)
    h2 = -pm(3) / 2. + pm(2) / 2.*1j
    h3 = -pm(3) / 2. - pm(2) / 2.*1j
    h = tridiag(h1, h2, h3, N)  # this creates 2N*2N hamiltonian that describes the system for a GIVEN k_3
    Mm = np.empty(shape=(0, 0))
    for ii in range(N):
        Mm=kronsum(Mm,pm(3)*(np.cos(b(ii,fs,N,rr)/2.)+2.))
        #Mm = kronsum(Mm, pm(3) * (b(ii, fs, N, rr) ))
    h=h+Mm
    w = LA.eigvalsh(h)  # this calculates the eigen values aka energies of the system for given k_3

    plt.plot(np.ones((np.size(w))) * k3[i]/np.pi, w, 'ko',
             markersize=0.3)  # this puts dots for each eigen value for a given k_3

fig2=plt.figure(2)


for ii in range(N):
    plt.plot(ii,b(ii,fs,N,rr)/np.pi,'ko',markersize=3)




fig3=plt.figure(3)


for ii in range(N):
    plt.plot(ii,(np.cos(b(ii,fs,N,rr)/2.)+2.+mg),'ko',markersize=3)

print(w)
plt.show()