import numpy as np  # this loads a library for linear algebra and comman math functions
import matplotlib.pyplot as plt  # this load library for ploting
from numpy import \
    linalg as LA  # this loads a library which is in the numpy for linearalgebra this is done for convience
from scipy import optimize as op


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


def ds(i, n, N,range):
    if n == 1:
        result = range
    elif n == 7:
        result = i
    elif n == 6:
        result = np.exp(-((i - N) / 50) ** 2)
    elif n == 2:
        if (i < N / 2):
            result = i/(N/2.)*range
        else:
            result = (N - i - 1.)/(N/2.)*range
    elif n==3:
        amount=N/2.5
        if (i < amount):
            result = i*range/(amount)
        elif (i >= amount)&(i < N-amount):
            result = range
        else :
            result=(N-i-1.)*range/(amount)
    elif n==4:
        amount=N/5
        if (i < amount):
            result = i*range/(amount)
        elif (i >= amount)&(i < N-amount):
            result = range
        else :
            result=(N-i-1.)*range/(amount)
    elif n==5:
        amount=N/10
        if (i < amount):
            result = i*range/(amount)
        elif (i >= amount)&(i < N-amount):
            result = range
        else :
            result=(N-i-1.)*range/(amount)

    return result
def deltaval(i,N,range):
    result=range
    return result
def sce(k1,lvl):

    k_step = 0.1

    k2 = 0.  # set k2 zero
    k3 = 0.  # and set k3 zero
    N = 150 # this is number of sites
    gamma = 0.  # m1 and m2 controls the distance between the weyl nodes
    m = -1.
    tx = 0.5
    t = 0.5
    k0 = np.pi/3.
    fs = 1
    mg = 0.9

    g0 = 1.
    T=0.
    idx=6





    # super conductivity



    h1 = -(m * (2. - np.cos(k3)) + 2 * tx * (np.cos(k1) - np.cos(k0) )) * pm(1) - 2 * t * np.sin(k3) * pm(2)
    h2 = +m * pm(1) / 2 + t * pm(3) / 1j
    h3 = +m * pm(1) / 2 - t * pm(3) / 1j
    H = tridiag(h1, h2, h3, N)  # this creates 2N*2N hamiltonian that describes the system for a GIVEN k_3

    w, U = LA.eigh(H)  # this calculates the eigen values aka energies of the system for given k_3
    # print(w[lvl])


    Pup=U[0:2*N-1:2,lvl]
    Pdown=U[1:2*N-1:2,lvl]
    # print(np.size(Pup))
    for i in range(4*N/10,6*N/10,1):
        plt.plot(i,abs(Pup[i]),'ko')
        plt.plot(i, abs(Pdown[i]), 'ro')
    fig1=plt.figure(1)
    # plt.ylim(0,np.max(abs(Pup)))



    return w[lvl]
print(sce(0.,150))
N=30
# for i in range(N):
#     plt.plot(i,np.real(sce(i,0)),'ko')
#     plt.plot(i, np.imag(sce(i, 0)), 'ro')
# root=op.newton(sce,-1.)
# print(root)
plt.show()