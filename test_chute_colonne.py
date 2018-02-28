#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:24:07 2018

@author: sdiop

Animation d'une colonne de sable qui tombe
"""

import main_verlet
import affichage

# paramètres :

# grains
R = 1
m = 1

# exterieur
L=20
e=2
y_min = -R
g = 1

# paramètres des contacts
Kc = 5e5/(2*R)
restitution = 0.3
Kd = restitution**2*Kc
Kd_2 = 0
Kc_p = 1e6/(2*R)
restitution_p = 0
Kd_p = restitution_p**2*Kc_p
Kt = 7
mu = 0.7
mu_p = 1
Kt_p = 7
alpha_0 = 0


dt = 1e-4

paroi = True
gravite = True




if __name__=='__main__':
    
    fig = plt.figure(1,figsize=(3.5,7))
    plt.axis("equal")
    ax = plt.axes(xlim=(-L/2,L/2), ylim=(0,2*L))
    fenetre = fig, ax
    grains = []
    for i in range(0,10):
        for j in range(0,5):
            grains.append(grain(R+(i-5)*2*R,R+j*2*R+(i%2)*R,0,0))
            
    ani = animer(grains)
#    X,Y,VX,VY,T = trajectoires(grains,0.1)
#    for k in range(len(X)):
#        plt.plot(X[k],Y[k],linewidth = 1)
    
    plt.show()
