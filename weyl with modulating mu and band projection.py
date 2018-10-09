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
def sce(mu,N,k0,mg,k_step,kstep2,wd,selection):

    if selection ==1:
        ws=0
        kstep2 = N * kstep2


        k1 = np.arange(-np.pi, np.pi, k_step)
        k2 = np.arange(-np.pi, np.pi, kstep2)
        k3 = np.arange(-np.pi, np.pi, k_step)


        mu_matrix = np.empty(shape=(0, 0))
        for ii in range(N):
            mu_matrix = kronsum(mu_matrix, pm(0) * mu[ii])
        kk=np.empty(shape=(3, 0))

        for i, val in enumerate(k1):
            for ii, val in enumerate(k2):
                for iii, val in enumerate(k3):



                    h1=-(np.cos(k3[iii])-2.+np.cos(k1[i])-np.cos(k0)-mg)*pm(1)-pm(3)*np.sin(k3[iii])
                    h3=-pm(2)/(2.*1j)-pm(1)/2.
                    h2=pm(2)/(2.*1j)-pm(1)/2.
                    H=tridiag(h1,h3,h2,N)

                    H=H-mu_matrix
                    hd=np.zeros((2*N,2*N))*1j
                    hd[0:2,2*N-2:2*N]=h2*np.exp(+1j*k2[ii])
                    hd[2*N-2:2*N,0:2]=h3*np.exp(-1j*k2[ii])
                    H=H+hd
                    # print(H)
                    # print('la')
                    # print(LA.eig(pm(0)))

                    w, U = LA.eigh(H)  # this calculates the eigen values aka energies of the system for given k_1
                    U = U[(w > -wd) & (w < wd), :]


                    w = w[(w > -wd) & (w < wd)]

                    if np.size(w)>0:
                        kk=np.append(kk,[[k1[i]],[k2[ii]],[k3[iii]]],1)
                        ws=ws+np.sum(w)
                        # plt.plot(np.ones((np.size(w))) * k1[i] ,np.ones((np.size(w))) * k2[ii] ,w, 'ko', markersize=0.8)


    elif selection ==2:
        k3 = 0.
        k2 = 0.

        # mu=-3.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
        # mu=mu+np.max(-mu)
        # mu=np.ones((N,))
        k1 = np.arange(-np.pi, np.pi, k_step)
        k2 = 0.
        k3 = 0.
        deltam = np.empty(shape=(0, 0))
        for ii in range(N):
            deltam = kronsum(deltam,pm(2)*1j*(-1*(np.cos(ds(ii,2,N,5)/2.)-np.cos(k0))))
        delta_up = deltam
        delta_down = delta_up.getH()
        delta_up = np.hstack((np.zeros((2 * N, 2 * N)),
                              delta_up))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
        delta_down = np.hstack((delta_down, np.zeros((2 * N, 2 * N))))
        delta_up = deltam
        delta_down = delta_up.getH()
        delta_up = np.hstack((np.zeros((2 * N, 2 * N)),
                              delta_up))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
        delta_down = np.hstack((delta_down, np.zeros((2 * N, 2 * N))))
        deltam = np.vstack((delta_up, delta_down))

        mu_matrix = np.empty(shape=(0, 0))
        for ii in range(N):
            mu_matrix = kronsum(mu_matrix, pm(0) * mu[ii])
        for i, val in enumerate(k1):
            h1 = -(np.cos(k3) - 2. + np.cos(k1[i]) - np.cos(k0) - mg) * pm(1) - pm(3) * np.sin(k3) * pm(3)
            h3 = -pm(2) / (2. * 1j) - pm(1) / 2.
            h2 = pm(2) / (2. * 1j) - pm(1) / 2.
            H = tridiag(h1, h3, h2, N)

            H = H - mu_matrix
            hd = np.zeros((2 * N, 2 * N)) * 1j
            hd[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(+1j * k2[i])
            hd[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(-1j * k2[i])
            H = H + hd
            h1m = -(np.cos(-k3) - 2. + np.cos(-k1[i]) - np.cos(k0) - mg) * pm(1) - pm(3) * np.sin(-k3) * pm(3)
            h3m = -pm(2) / (2. * 1j) - pm(1) / 2.
            h2m = pm(2) / (2. * 1j) - pm(1) / 2.
            Hm = tridiag(h1m, h3m, h2m, N)
            hdm = np.zeros((2 * N, 2 * N)) * 1j
            hdm[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(-1j * k2[i])
            hdm[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(+1j * k2[i])
            Hm = Hm - mu_matrix + hdm
            h_up = np.hstack((H, np.zeros((2 * N,
                                           2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
            h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
            h_bog = np.vstack((h_up, h_down))

            h_bog = deltam + h_bog
            w, U = LA.eigh(h_bog)
            U = U[(w > -1) & (w < 1), :]
            w = w[(w > -1) & (w < 1)]
            ax.plot(np.ones((np.size(w))) * k2[i] / np.pi, w, 'ko', markersize=0.3)

    return kk,ws

N=10
k0 =np.pi/3.
mg = 0.9*0.
k_step = 0.1
kstep2=0.1
mu=-2.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
mu=mu+np.max(-mu)
mu=np.ones((N,))*1.5
wd=0.2
# fig1 = plt.figure(1)
# ax = fig1.gca(projection = '3d')
print(mu)
k11,ws=sce(mu,N,k0,mg,k_step,kstep2,wd,1)
NN=np.size(k11)
print(np.size(k11))
print(ws)
# k22=np.array(-N/4,N/4,1)
# # k22=k_step*k22
# print(k11)
# print(np.diff(k11))
plt.show()