import numpy as np  # this loads a library for linear algebra and comman math functions
import matplotlib.pyplot as plt  # this load library for ploting
from numpy import \
    linalg as LA  # this loads a library which is in the numpy for linearalgebra this is done for convience
from scipy import optimize as op
# import multiprocessing as mp


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
def eigenvecfunc(k1,k2,k3,k0,mg,deltam,mu_matrix):



    h1 = -(np.cos(k3) - 2. + np.cos(k1) - np.cos(k0) - mg) * pm(1) - pm(3) * np.sin(k3) * pm(3)
    h3 = -pm(2) / (2. * 1j) - pm(1) / 2.
    h2 = pm(2) / (2. * 1j) - pm(1) / 2.
    H = tridiag(h1, h3, h2, N)

    H = H - mu_matrix
    hd = np.zeros((2 * N, 2 * N)) * 1j
    hd[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(+1j * k2)
    hd[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(-1j * k2)
    H = H + hd
    h1m = -(np.cos(-k3) - 2. + np.cos(-k1) - np.cos(k0) - mg) * pm(1) - pm(3) * np.sin(-k3) * pm(
        3)
    h3m = -pm(2) / (2. * 1j) - pm(1) / 2.
    h2m = pm(2) / (2. * 1j) - pm(1) / 2.
    Hm = tridiag(h1m, h3m, h2m, N)
    hdm = np.zeros((2 * N, 2 * N)) * 1j
    hdm[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(-1j * k2)
    hdm[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(+1j * k2)
    Hm = Hm - mu_matrix + hdm
    h_up = np.hstack((H, np.zeros((2 * N,
                                   2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
    h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
    h_bog = np.vstack((h_up, h_down))

    h_bog = deltam + h_bog
    w, U = LA.eigh(h_bog)
    return w


def kcutoff(mu, N, k0, mg, k_step, wd, selection):
    if selection == 1:
        fig2 = plt.figure(2)
        k3 = 0.
        k2 = 0.

        # mu=-3.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
        # mu=mu+np.max(-mu)
        # mu=np.ones((N,))
        k2 = np.arange(-np.pi, np.pi, k_step)
        k1 = np.ones((np.size(k2),)) * k0
        k1 = np.arange(-np.pi, np.pi, k_step)
        k2 = np.ones((np.size(k1),)) * 0.
        k3 = 0.

        mu_matrix = np.empty(shape=(0, 0))
        for ii in range(N):
            mu_matrix = kronsum(mu_matrix, pm(0) * mu[ii])
        k11 = np.zeros((np.size(k1),))

        for i, val in enumerate(k1):

            h1 = -(np.cos(k3) - 2. + np.cos(k1[i]) - np.cos(k0) - mg) * pm(1) - pm(3) * np.sin(k3) * pm(3)
            h3 = -pm(2) / (2. * 1j) - pm(1) / 2.
            h2 = pm(2) / (2. * 1j) - pm(1) / 2.
            H = tridiag(h1, h3, h2, N)

            H = H - mu_matrix
            hd = np.zeros((2 * N, 2 * N)) * 1j
            hd[0:2, 2 * N - 2:2 * N] = h2 * np.exp(+1j * k2[i])
            hd[2 * N - 2:2 * N, 0:2] = h3 * np.exp(-1j * k2[i])
            H = H + hd
            # print(H)
            # print('la')
            # print(LA.eig(pm(0)))

            w, U = LA.eigh(H)  # this calculates the eigen values aka energies of the system for given k_1
            U = U[(w > -wd) & (w < wd), :]

            w = w[(w > -wd) & (w < wd)]
            if np.size(w) > 0:
                k11[i] = k1[i]

            plt.plot(np.ones((np.size(w))) * k1[i] / np.pi, w, 'ko', markersize=0.8)


    elif selection == 2:
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
            deltam = kronsum(deltam, pm(2) * 1j * (-1 * (np.cos(ds(ii, 2, N, 5) / 2.) - np.cos(k0))))
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
            # plt.plot(np.ones((np.size(w))) * k2[i] / np.pi, w, 'ko', markersize=0.3)

    k11 = k11[k11 != 0]

    return k11
def sce(delta,N,mu,k1,k2,k3,k_step):
    # pool=mp.Pool(processes=8)
    # N=10
    # k3=0.
    # k2=0.
    g0=-10.
    k0=np.pi/3.
    mg=0.9*0.

    #mu=-10.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
    #mu=mu+np.max(-mu)
    #mu=5.*np.ones((N,))
    # k1 = np.arange(-np.pi, np.pi, k_step)
    # k2 = np.arange(-np.pi, np.pi, k_step)
    # k3 = np.arange(-np.pi, np.pi, k_step)
    V=np.size(k1)*np.size(k2)*np.size(k3)
    V=np.size(np.arange(-np.pi, np.pi, k_step))**3

    smt = 0.
    deltam = np.empty(shape=(0, 0))
    for ii in range(N):
        deltam = kronsum(deltam, pm(2) * 1j * delta[ii])
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
        for jj, val in enumerate(k2):
            for ll, val in enumerate(k3):

                h1=-(np.cos(k3[jj])-2.+np.cos(k1[i])-np.cos(k0)-mg)*pm(1)-pm(3)*np.sin(k3[ll])*pm(3)
                h3=-pm(2)/(2.*1j)-pm(1)/2.
                h2=pm(2)/(2.*1j)-pm(1)/2.
                H=tridiag(h1,h3,h2,N)

                H=H-mu_matrix
                hd=np.zeros((2*N,2*N))*1j
                hd[0:2,int(2*N-2):2*N]=h2*np.exp(+1j*k2[jj])
                hd[int(2*N-2):2*N,0:2]=h3*np.exp(-1j*k2[jj])
                H=H+hd
                h1m=-(np.cos(-k3[ll])-2.+np.cos(-k1[i])-np.cos(k0)-mg)*pm(1)-pm(3)*np.sin(-k3[ll])*pm(3)
                h3m=-pm(2)/(2.*1j)-pm(1)/2.
                h2m=pm(2)/(2.*1j)-pm(1)/2.
                Hm=tridiag(h1m,h3m,h2m,N)
                hdm = np.zeros((2 * N, 2 * N)) * 1j
                hdm[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(-1j * k2[jj])
                hdm[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(+1j * k2[jj])
                Hm=Hm-mu_matrix+hdm
                h_up = np.hstack((H, np.zeros((2 * N,
                                               2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
                h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
                h_bog = np.vstack((h_up, h_down))

                h_bog = deltam + h_bog
                w, U = LA.eigh(h_bog)  # this calculates the eigen values aka energies of the system for given k_3
                # U=U[(w>-wd)&(w<wd),:]
                #plt.plot(np.ones((np.size(w))) * k1[i] / np.pi, w, 'ko', markersize=0.3)
                n = np.arange(0, 2 * N - 1, 1)
                #print(w)
                #print('la')
                U = np.fliplr(U)
                u=U[0,:]
                # ph=np.angle(u)
                #
                # ph=np.squeeze(np.asarray(ph))
                #print(np.diag(ph))
                # w = np.flipud(w)
                # idx=np.arange(0,N-1,1)
                # u2k=U[1:2*N:2, 0:2*N]
                #
                # cv1k = np.conj(U[(2 * N):(4*N-1):2, 0:2*N])
                # #u2k=U[1:2*N:2, :]
                #
                # #cv1k = np.conj(U[(2 * N):(4*N-1):2, :])
                # #U=U
                # #U=U*np.exp(1j*np.diag(np.squeeze(np.asarray(np.angle(np.multiply(u2k[0,:],cv1k[0,:]))))))
                # #print(np.size(np.exp(-1j*np.diag(np.squeeze(np.asarray(np.angle(np.multiply(u2k[0,:],cv1k[0,:]))))/2))))
                # #print('la')
                u2k=U[1:2 * N:2, 0:2 * N]
                cv1k = np.conj(U[(2 * N):(4 * N - 1):2, 0:2 * N])

                # aaa=u2k
                # print(np.angle(np.multiply(u2k[0,:],cv1k[0,:]))/np.pi)

                smt=np.sum(np.multiply(u2k,cv1k),axis=1)+smt
                #print(smt)

                #print('la')

                # for idx in range(N):
                #     for j in n:
                #         un2ki = U[2 * idx - 1, n[j]]
                #         cv1nki = np.conj(U[2 * N + 2 * (idx - 1), n[j]])
                #         #En = w[n[j]]
                #         smt[idx] = smt[idx] + un2ki * cv1nki
                #         #smt=np.real(smt)
                #     smt[idx] = smt[idx] * g0 * k_step**3 / (2. * np.pi)**3 - 1 * delta[idx]
                    #smt[idx] = smt[idx]  - 1 * delta[idx]
    #print(smt* g0 * k_step**3 / (2. * np.pi)**3 - 1 * delta)
    smt=np.squeeze(np.asarray(smt))
    return np.real(smt)* g0/V - 1. * delta

N=2
k0 =np.pi/3.
mg = 0.9*0.
wd=1
delta=5.*np.ones((N,))
mu=-3.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
mu=np.zeros(N,)
for i in range(N):
    mu[i]=ds(i,2,N,10.)
mu=-mu+np.max(mu)
# mu=3.*np.ones((N,))
mu=np.ones((N,))*2.5
k_step=0.5
k1 = kcutoff(mu, N, k0, mg, k_step, wd, 1)
NN=np.size(k1)
k2=np.arange(-NN/2.,NN/2.,1)
k2=k_step*k2
k3 = k2
print(np.size(k1)*np.size(k2)*np.size(k3))



# root=op.fsolve(sce,delta,args=(N,mu,k1,k2,k3,k_step))
# print(root)
# print('yukardaki root')
# print(10.*np.sin(np.pi*np.arange(0,N,1)/N))
# print('yukardaki mu')
# print(sce(delta,N,mu,k1,k2,k3,k_step))
# print('yukardaki sce')
# k=np.arange(0,N,1)
# fig1=plt.figure(1)
# plt.plot(k,mu,'ko-',k,root,'ro-')
# plt.ylim([0,20])
n=np.arange(0.,3.,0.1)
rr=n*0.
for jj in range(np.size(n)):
    mu = np.ones((N, ))*n[jj]
    k_step = 0.8
    k1 = kcutoff(mu, N, k0, mg, k_step, wd, 1)
    NN = np.size(k1)
    k2 = np.arange(-NN / 2., NN / 2., 1)
    k2 = k_step * k2
    k3 = k2
    root = op.fsolve(sce, delta, args=(N, mu, k1, k2, k3, k_step))
    rr[jj]=np.max(root)
    print(np.max(root))
fig4=plt.figure(4)
plt.plot(n,rr,'ko-')
plt.ylim((0, np.max(mu)))
plt.ylabel(r'$|\Delta(\mu)|$')
plt.xlabel(r'$\mu$')
plt.xticks([0, 0.5, 1, 1.5, 2,2.5,3],
               [r'$0$', r'$0.5$', r'$1$', r'$1.5$', r'$2$',r'$2.5$',r'$3$'])
plt.yticks([0, 1, 2, 3,4,5],
               [r'$0$', r'$1$', r'$2$', r'$3$', r'$4$',r'$5$'])

print(sce(root,N,mu,k1,k2,k3,k_step))
print('yukardaki sce with root')

plt.show()
# print(sce(delta,N))
# print('yukardaki sce')
# print(delta)
# print('yukardaki delta')