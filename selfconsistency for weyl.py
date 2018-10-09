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
def sce(delta):

    k_step = 0.1
    k1 = np.arange(-np.pi, np.pi, k_step)  # this kreates k1 momentum array from -pi to pi
    k2 = 0.  # set k2 zero
    k3 = 0.  # and set k3 zero
    N = 30  # this is number of sites
    gamma = 0.  # m1 and m2 controls the distance between the weyl nodes
    m = -1.
    tx = 0.5
    t = 0.5
    k0 = 0.
    fs = 1
    mg = 0.9

    g0 = 1.
    T=0.
    idx=6





    # super conductivity
    smt = 0.
    for i, val in enumerate(k1):
        if T!=0. :
            h1 = -(m * (2. - np.cos(k3)) + 2 * tx * (np.cos(k1[i]) - np.cos(k0) - mg)) * pm(1) - 2 * t * np.sin(k3) * pm(2)
            h2 = +m * pm(1) / 2 + t * pm(3) / 1j
            h3 = +m * pm(1) / 2 - t * pm(3) / 1j
            H = tridiag(h1, h2, h3, N)  # this creates 2N*2N hamiltonian that describes the system for a GIVEN k_3
            hm1 = -(m * (2. - np.cos(-k3)) + 2 * tx * (np.cos(-k1[i]) - np.cos(k0) - mg)) * pm(1) - 2 * t * np.sin(
                -k3) * pm(2)
            hm2 = +m * pm(1) / 2 + t * pm(3) / 1j
            hm3 = +m * pm(1) / 2 - t * pm(3) / 1j
            Hm = tridiag(hm1, hm2, hm3, N)
            h_up = np.hstack((H, np.zeros((2 * N,
                                           2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
            h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
            h_bog = np.vstack((h_up, h_down))

            delta = np.empty(shape=(0, 0))

            for ii in range(N):
                delta = kronsum(delta, pm(2) * 1j * (-1 * (np.cos(ds(ii, fs, N, rr) / 2.) - np.cos(k0))))
            # delta_up=np.kron(np.identity(N),delta_0*pm(2)*1j/slope)
            delta = np.empty(shape=(0, 0))
            for ii in range(N):
                delta = kronsum(delta, pm(2) * 1j * deltaval(ii, N, rr))
            # delta_up=np.kron(np.identity(N),delta_0*pm(2)*1j/slope)
            delta_up = delta
            delta_down = delta_up.getH()
            delta_up = np.hstack((np.zeros((2 * N, 2 * N)),
                                  delta_up))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
            delta_down = np.hstack((delta_down, np.zeros((2 * N, 2 * N))))
            delta = np.vstack((delta_up, delta_down))
            h_bog = delta + h_bog
            w, U = LA.eigh(h_bog)  # this calculates the eigen values aka energies of the system for given k_3
            n = np.arange(0, 2 * N - 1, 1)

            for j in n:
                un1ki = U[2 * idx - 2, n[j]]
                un2ki = U[2 * idx - 1, n[j]]
                cv2nki = np.conj(U[2 * N + 2 * (idx - 1) + 1, n[j]])
                cv1nki = np.conj(U[2 * N + 2 * (idx - 1), n[j]])
                En = w[n[j]]
                smt = smt + cv2nki * un1ki * (1. + np.exp(En / T) ** (-1.)) + un2ki * cv1nki * (
                            1. - (1. + np.exp(En / T) ** (-1.)))
        elif T==0:
            h1 = -(m * (2. - np.cos(k3)) + 2 * tx * (np.cos(k1[i]) - np.cos(k0) - mg)) * pm(1) - 2 * t * np.sin(
                k3) * pm(2)
            h2 = +m * pm(1) / 2 + t * pm(3) / 1j
            h3 = +m * pm(1) / 2 - t * pm(3) / 1j
            H = tridiag(h1, h2, h3, N)  # this creates 2N*2N hamiltonian that describes the system for a GIVEN k_3
            hm1 = -(m * (2. - np.cos(-k3)) + 2 * tx * (np.cos(-k1[i]) - np.cos(k0) - mg)) * pm(1) - 2 * t * np.sin(
                -k3) * pm(2)
            hm2 = +m * pm(1) / 2 + t * pm(3) / 1j
            hm3 = +m * pm(1) / 2 - t * pm(3) / 1j
            Hm = tridiag(hm1, hm2, hm3, N)
            h_up = np.hstack((H, np.zeros((2 * N,
                                           2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
            h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
            h_bog = np.vstack((h_up, h_down))

            delta = np.empty(shape=(0, 0))

            for ii in range(N):
                delta = kronsum(delta, pm(2) * 1j * (-1 * (np.cos(ds(ii, fs, N, rr) / 2.) - np.cos(k0))))
            # delta_up=np.kron(np.identity(N),delta_0*pm(2)*1j/slope)
            delta = np.empty(shape=(0, 0))

            for ii in range(N):
                delta = kronsum(delta, pm(2) * 1j * deltaval(ii,N,rr))
            # delta_up=np.kron(np.identity(N),delta_0*pm(2)*1j/slope)

            delta_up = delta
            delta_down = delta_up.getH()
            delta_up = np.hstack((np.zeros((2 * N, 2 * N)),
                                  delta_up))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
            delta_down = np.hstack((delta_down, np.zeros((2 * N, 2 * N))))
            delta = np.vstack((delta_up, delta_down))
            h_bog = delta + h_bog
            w, U = LA.eigh(h_bog)  # this calculates the eigen values aka energies of the system for given k_3
            n = np.arange(0, 2 * N - 1, 1)
            U=np.fliplr(U)
            w=np.flipud(w)

            for j in n:
                un1ki = U[2 * idx - 2, n[j]]
                un2ki = U[2 * idx - 1, n[j]]
                cv2nki = np.conj(U[2 * N + 2 * (idx - 1) + 1, n[j]])
                cv1nki = np.conj(U[2 * N + 2 * (idx - 1), n[j]])
                En = w[n[j]]
                smt = smt  + un2ki * cv1nki




    smt = smt * g0 * k_step / (2. * np.pi)-1*deltaval(idx,N,rr)
    return smt
print(sce(0.05))
N=30
# for i in range(N):
#     plt.plot(i,np.real(sce(i,0)),'ko')
#     plt.plot(i, np.imag(sce(i, 0)), 'ro')
# root=op.newton(sce,-1.)
# print(root)
# plt.show()