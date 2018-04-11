#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:13:06 2018

@author: sdiop
"""

# paramètres :

# exterieur
R = 1
m = 1
L=200000
e=20000
y_min = -900000

# paramètres des contacts
Kc = 1e4
rest = 0.8
Kd = rest**2 * Kc

Kc_p = 1e4
rest_p = 0.4
Kd_p = rest_p**2 * Kc_p


Kt = 0
mu = 0.7
Kt_p = 0
mu_p = 0.7


dt = 1e-3

paroi = False
gravite = False
remove= False




if __name__=='__main__':
    m = 1
    
    #%%
    
    eps = 0
    for k in range(1,40):
        grains = []
        v0 = 1e-6* (k+3)**4
        grains.append(grain(-1.001,eps-10,v0,0))
        grains.append(grain(1,-10,0,0))
        
        T,X,Y,VX,VY = trajectoires(grains, 5)
        
        if (X[1][-1]-X[0][-1])**2+(Y[1][-1]-Y[0][-1])**2 < 4 : col = 'red'
        else : col = 'blue'
        E = (VX[0]**2 + VY[0]**2 + VX[1]**2 + VY[1]**2)/2 
        E = E/E[0]
        perte = 1-E[-1]
        plt.scatter(VX[0][0], perte, color = col)   
        
    plt.title(r"Perte relative d'$E_{c}$ en fonction de la vitesse initiale")
    plt.xlabel(r'$v_{0}$')
    plt.ylabel(r'$\eta$')
    plt.grid()
    plt.show()
    
    
    #%%
    eps = 0
    v0 = 1
    y0 = 1e6
    x0 = 1e6
    grains = [grain(-1.1+x0,eps+y0,v0,0), grain(1+x0,0+y0,0,0)]
    
    fig = plt.figure(1,figsize=(6,6))
    plt.axis("equal")
    ax = plt.axes(xlim=(-5+x0,5+x0), ylim=(-5+y0,5+y0))
    fenetre = fig, ax
    ani = animer(grains)
    plt.show()
    
    #%%
    eps = 0
    v0 = 1
    
    for n in range(10):
        y0 = 10**(n-7)
        x0 = 10**(n-7)
        grains = [grain(-1.1+x0,eps+y0,v0,0), grain(1+x0,0+y0,0,0)]
        T,X,Y,VX,VY = trajectoires(grains, 5)
        E = (VX[0]**2 + VY[0]**2 + VX[1]**2 + VY[1]**2)/2 
        E = E/E[0]
    #    plt.plot(E)
    #    plt.show()
        plt.scatter(y0,VY[1][-1])
        plt.xscale('log')
    plt.show()
    