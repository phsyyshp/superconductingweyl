import numpy as np
import matplotlib  as plt
def cn(kx,ky,n):
    tx=1
    ty=1
    m=1
    k0=np.pi/3
    A=-2*ty*np.sin(ky)
    B=2*tx*(np.cos(kx)-np.cos(k0))+m*(1-np.cos(ky))
    E1=np.sqrt(B**2+A**2)
    if n==3:

        if np.sqrt(2*E1*(E1-B))==0:
            result=0
        else :
            result = 1j * (-B + np.sqrt(B ** 2 + A ** 2)) / np.sqrt(2 * E1 * (E1 - B))
    elif n==1:

        if np.sqrt(2*E1*(E1+B))==0:
            result=0
        else :
            result = -1j * (B + np.sqrt(B ** 2 + A ** 2)) / np.sqrt(2 * E1 * (E1 + B))
    elif n==4:

        if np.sqrt(2*E1*(E1-B))==0:
            result=0
        else :
            result = A / np.sqrt(2 * E1 * (E1 - B))
    elif n==2:

        if np.sqrt(2*E1*(E1+B))==0:
            result=0
        else :
            result = A / np.sqrt(2 * E1 * (E1 + B))
    return result

def s(qx,qy,qxp,qyp):
    inc=0.4
    kx=np.arange(-np.pi,np.pi,inc)
    ky = np.arange(-np.pi, np.pi, inc)
    kxb = np.arange(-np.pi, np.pi, inc)
    kyb = np.arange(-np.pi, np.pi, inc)
    result=0.
    for ikx in kx:
        for iky in ky:
            for ikyb in kyb:
                for ikxb in kxb:
                    result=result+(np.conj(cn(qx+ikx,qy+iky,3))*cn(ikx,iky,1)-np.conj(cn(qx+ikx,qy+iky,4))*cn(ikx,iky,2))*(np.conj(cn(qxp+ikxb,qyp+ikyb,1))*cn(ikxb,ikyb,3)-np.conj(cn(qxp+ikxb,qyp+ikyb,2))*cn(ikxb,ikyb,4))

    return result
print(s(qx,qy,qxp,qyp))