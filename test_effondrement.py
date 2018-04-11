#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:24:07 2018

@author: sdiop

Animation d'une colonne de sable qui tombe
"""


# paramètres :

# exterieur
R = 1
m = 1
L=20000000
e=0
g = 1
y_min = -200000

# paramètres des contacts
Kc = 1e4
rest = 0.8
Kd = rest**2 * Kc

Kc_p = 1e4
rest_p = 0.4
Kd_p = rest_p**2 * Kc_p
Kt = 0
mu = 0.7
mu_p = 0.7
Kt_p = 10


dt = 1e-3

paroi = True
gravite = True
remove = False

if __name__=='__main__':
    m=1

    
    #%%
    
    grains = []
    for i in range(1,10):
        for j in range(0,5):
            grains.append(grain(R+(i-5)*2*R,R+j*2*(R+0.2)+(i%2)*R,1*(rd.random()-0.5),1*(rd.random()-1)))
    fig = plt.figure(1,figsize=(10,10))
    plt.axis("equal")
    ax = plt.axes(xlim=(-30,30), ylim=(0,60))
    fenetre = fig, ax
    ani = animer(grains)
    plt.grid()
    plt.show()
    
    
    #%%
    
    grains = []
    for i in range(1,10):
        for j in range(0,5):
            grains.append(grain(R+(i-5)*2*R,R+j*2*(R+0.2)+(i%2)*R,1*(rd.random()-0.5),1*(rd.random()-1)))

    temps = 30
    T,X,Y,VX,VY = trajectoires(grains, temps)
    Ec = 0
    Ep = 0
    for k in range(np.shape(X)[0]):
        Ec += (VX[k]**2 + VY[k]**2)/2
        Ep += Y[k]
    E = Ec+Ep
    for k in range(len(grains)):
        plt.subplot(1,2,1)
        plt.plot(X[k],Y[k])
    plt.subplot(1,2,2)
    plt.plot(T,E/E[0],label = 'E')
    plt.plot(T,Ec/E[0],label = 'E_c')
    plt.plot(T,Ep/E[0],label = 'E_p')
    plt.xlim(0,temps)
    plt.legend()
    plt.grid()
    plt.show()
    

    