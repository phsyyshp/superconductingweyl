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


def kcutoff(mu,N,k0,mg,k_step,kstep2,wd,selection):
    if selection ==1:
        ws = 0.
        kstep2=N*kstep2


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



                    h1=-(np.cos(k3[iii])-2.+np.cos(k1[i])-np.cos(k0)-mg)*pm(1)-np.sin(k3[iii])*pm(3)
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
                        # ws=ws+nw
                    # plt.plot(np.ones((np.size(w))) * k1[i] , w, 'ko', markersize=0.8)


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
            plt.plot(np.ones((np.size(w))) * k2[i] / np.pi, w, 'ko', markersize=0.3)


    return kk , ws

def sce(delta,N,mu,kk,k_step,kstep2,wdo):
    # pool=mp.Pool(processes=8)
    # N=10
    # k3=0.
    # k2=0.
    g0=-100.
    k0=np.pi/3.
    mg=0.9*0.
    kstep2=N*kstep2

    #mu=-10.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
    #mu=mu+np.max(-mu)
    #mu=5.*np.ones((N,))
    # k1 = np.arange(-np.pi, np.pi, k_step)
    # k2 = np.arange(-np.pi, np.pi, k_step)
    # k3 = np.arange(-np.pi, np.pi, k_step)

    V=np.size(np.arange(-np.pi, np.pi, k_step))**2*np.size(np.arange(-np.pi, np.pi, kstep2))

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
    for i in range(np.size(kk,1)):


        h1=-(np.cos(kk[2,i])-2.+np.cos(kk[0,i])-np.cos(k0)-mg)*pm(1)-pm(3)*np.sin(kk[2,i])
        h3=-pm(2)/(2.*1j)-pm(1)/2.
        h2=pm(2)/(2.*1j)-pm(1)/2.
        H=tridiag(h1,h3,h2,N)

        H=H-mu_matrix
        hd=np.zeros((2*N,2*N))*1j
        hd[0:2,int(2*N-2):2*N]=h2*np.exp(+1j*kk[1,i])
        hd[int(2*N-2):2*N,0:2]=h3*np.exp(-1j*kk[1,i])
        H=H+hd
        h1m=-(np.cos(-kk[2,i])-2.+np.cos(-kk[0,i])-np.cos(k0)-mg)*pm(1)-pm(3)*np.sin(-kk[2,i])
        h3m=-pm(2)/(2.*1j)-pm(1)/2.
        h2m=pm(2)/(2.*1j)-pm(1)/2.
        Hm=tridiag(h1m,h3m,h2m,N)
        hdm = np.zeros((2 * N, 2 * N)) * 1j
        hdm[0:2, int(2 * N - 2):2 * N] = h2 * np.exp(-1j * kk[1,i])
        hdm[int(2 * N - 2):2 * N, 0:2] = h3 * np.exp(+1j * kk[1,i])
        Hm=Hm-mu_matrix+hdm
        h_up = np.hstack((H, np.zeros((2 * N,
                                       2 * N))))  # this following three lines creates a block diagonal matrix one block is H and the other block is -H^t
        h_down = np.hstack((np.zeros((2 * N, 2 * N)), -np.transpose(Hm)))
        h_bog = np.vstack((h_up, h_down))

        h_bog = deltam + h_bog
        w, U = LA.eigh(h_bog)  # this calculates the eigen values aka energies of the system for given k_3

        n = np.arange(0, 2 * N - 1, 1)
        wd = np.sqrt(wdo**2+(np.sum(delta)/np.size(delta))**2)
        # wd = np.sqrt(wdo ** 2 + (np.min(abs(delta))) ** 2)
        # wd = np.sqrt(wdo ** 2 + 5 ** 2)
        # print(np.sum(delta)/np.size(delta))

        #

        # U = U[:,(w > -wd) & (w < wd)]
        # print(U.shape)
        # w = w[(w > -wd) & (w < wd)]

        # U = np.fliplr(U)

        #mm=np.size(w)
        # print(mm)
        # mm=mm/4


        # #print('la')
        # print(mm)
        # u2k=U[1:2 * N:2, 0:2 * mm]
        # cv1k = np.conj(U[(2 * N):(4 * N - 1):2, 0:2 * mm])

        # smt=np.sum(np.multiply(u2k,cv1k),axis=1)+smt



        if np.size(w[(w > 0) & (w < wd)])>0:
            # print(w)
            # print(w[(w > -wd) & (w < wd)])
            # print(np.size(w[(w > -wd) & (w < wd)]))
            U = U[:,(w > 0) & (w < wd)]
            # print(U.shape)
            w = w[(w > 0) & (w < wd)]

            # U = np.fliplr(U)

            mm=np.size(w)
            # mm=mm/4
            # print(mm)


            # #print('la')
            # print(mm)
            u2k=U[2:2 * N-1:2, :]
            cv1k = np.conj(U[(2 * N):(4 * N - 1):2, :])

            smt=np.sum(np.multiply(u2k,cv1k),axis=1)+smt

    smt=np.squeeze(np.asarray(smt))
    print(smt*g0/V)
    return np.real(smt)* g0/V - 1. * delta

N=2
k0 =np.pi/3.
mg = 0.9*0.
wd=0.2
delta=3.*np.ones((N,))
mu=-3.*np.sin(np.pi*np.arange(0,N,1)/(N-1))
mu=np.zeros(N,)
for i in range(N):
    mu[i]=ds(i,2,N,0.3)
# mu=-mu+np.max(mu)
# mu=0.3*np.ones((N,))
mu=np.ones((N,))*0.3
k_step=0.2
kstep2=0.1
print(mu)
kk, ws = kcutoff(mu, N, k0, mg, k_step,kstep2 ,wd, 1)
sizeofks=np.size(kk)/3
print(sizeofks)

print(ws)
n=np.arange(0.01,0.3,0.01)
rr=n*0.
g0=-350.
# for jj in range(np.size(n)):
#     mu = np.ones((N, ))*n[jj]
#     k_step = 0.2
#     kstep2 = 0.1
#     kk, ws = kcutoff(mu, N, k0, mg, k_step, kstep2, wd, 1)
#
#     root = op.fsolve(sce, delta, args=(N,mu,kk,k_step,kstep2,wd))
#     rr[jj]=np.max(root)
#     print(np.max(root))
# fig4=plt.figure(4)
# plt.plot(n,rr,'ko-')
# wd=0.2
# vf=np.sin(np.pi/3)
# delta=5*wd*np.exp(-6*2*np.pi**2*vf/n**2/(-g0))
# plt.plot(n,delta,'ro-')
# plt.ylim((0, np.max(rr)))
# plt.ylabel(r'$|\Delta(\mu)|$',fontsize=25)
# plt.xlabel(r'$\mu$',fontsize=25)
# print(ws/sizeofks)
#
root=op.fsolve(sce,delta,args=(N,mu,kk,k_step,kstep2,wd))
print(root)
print('yukardaki root')
print(mu)
print('yukardaki mu')
print(sce(root,N,mu,kk,k_step,kstep2,wd))
print('yukardaki sce')
k=np.arange(0,N,1)
fig1=plt.figure(1)
plt.plot(k,mu,'ko-',label=r'$\mu(\alpha_i)$')
fig2=plt.figure(2)
print(root)
print('bu root')
plt.plot(k,root,'ro-')
plt.xlabel(r'$\frac{\alpha_i}{n}$',fontsize=40)
plt.ylim([np.min(root)-1,np.max(root)+1])
# # plt.xticks([0, 2.5, 5,7, 9],[r'$0$', r'$2.5$', r'$5$', r'$7$', r'$9$'],fontsize=25)
# # plt.yticks([0, 2, 4,6, 8,10],[r'$0$', r'$2$', r'$4$', r'$6$', r'$8$', r'$10$'],fontsize=25)
# # plt.ylim([0,20])
# # n=np.arange(0.,2.,0.05)
# # rr=n*0.
# # for jj in range(np.size(n)):
# #     mu = np.ones((N, ))*n[jj]
# #     k_step = 0.3
# #     kstep2 = 0.15
# #     kk, ws = kcutoff(mu, N, k0, mg, k_step, kstep2, wd, 1)
# #
# #     root = op.fsolve(sce, delta, args=(N,mu,kk,k_step,kstep2))
# #     rr[jj]=np.max(root)
# #     print(np.max(root))
# # fig4=plt.figure(4)
# # plt.plot(n,rr,'ko-')
# # plt.ylim((0, np.max(rr)))
# # plt.ylabel(r'$|\Delta(\mu)|$',fontsize=25)
# # plt.xlabel(r'$\mu$',fontsize=25)
# # plt.xticks([0, 1, 2, 3,4,5],
# #                [r'$0$', r'$1$', r'$2$', r'$3$', r'$4$',r'$5$'],fontsize=25)
# # plt.yticks([0, 1, 2, 3,4,5],
# #                [r'$0$', r'$1$', r'$2$', r'$3$', r'$4$',r'$5$'],fontsize=25)
# #
# # print(sce(root,N,mu,kk,k_step))
# # print('yukardaki sce with root')
#
#
#
# ## wd plot
# # n=np.arange(0.,3.,0.1)
# # rr=n*0.
# # mu=1.*np.ones((N,))
# # for jj in range(np.size(n)):
# #     wd = n[jj]
# #     k_step = 0.3
# #     kstep2 = 0.15
# #     kk, ws = kcutoff(mu, N, k0, mg, k_step, kstep2, wd, 1)
# #
# #     root = op.fsolve(sce, delta, args=(N,mu,kk,k_step,kstep2,wd))
# #     rr[jj]=np.max(root)
# #     print(np.max(root))
# # fig4=plt.figure(4)
# # plt.plot(n,rr,'ko-')
# # plt.ylim((0, np.max(rr)))
# # plt.ylabel(r'$|\Delta(\omega)|$',fontsize=25)
# # plt.xlabel(r'$\omega$',fontsize=25)
plt.show()
# print(sce(delta,N))
# print('yukardaki sce')
# print(delta)
# print('yukardaki delta')