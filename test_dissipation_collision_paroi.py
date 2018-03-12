#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 20:41:45 2018

@author: sdiop
"""

# grains
R = 1
m = 1
L=40
e=0
y_min = -R
g = 1

# paramètres des contacts
Kc = 1e4 * m*g/R
Kt = 7
mu = 0.5
gamma = 200
xi = 10

Kc_p = 1e4*m*g/R
Kt_p = 7
mu_p = 0.5
gamma_p = 20
xi_p = 10


dt = 1e-3
paroi = True
gravite = False
remove = False


if __name__=='__main__':
    
    grains = [grain(3*k+0.1,2,0.,-(k+4)) for k in range(-3,4)]
    
    #%%
    
    T,X,Y,VX,VY,FX,FY,FpX,FpY,FvX,FvY,G = trajectoires(grains, 1)
    for k in range(len(grains)):
        plt.scatter(np.abs(VY[k][0]),1-VY[k][-1]**2/VY[k][0]**2)  
        
    plt.title("Perte relative d'Ec en fonction de la vitesse initiale")
    

    plt.legend()
    plt.grid()
    plt.show()
    
    
    #%%
    grains = [grain(3*k+0.1,2,0,-(k+4)) for k in range(-3,4)]
    
    fig = plt.figure(1,figsize=(6,6))
    plt.axis("equal")
    ax = plt.axes(xlim=(-11,11), ylim=(0,22))
    fenetre = fig, ax
    ani = animer(grains)
    plt.show()