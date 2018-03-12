#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:13:06 2018

@author: sdiop
"""

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
Kc = 1e4 * m*g/R
restitution = 0.3
Kd = restitution**2*Kc
Kc_p = 1e4 * m*g/R
restitution_p = 0.3
Kd_p = restitution_p**2*Kc_p
Kt = 10
mu = 0.7
mu_p = 0.5
Kt_p = 7

gamma = 300
xi = 10


dt = 1e-4
paroi = False
gravite = False
remove=False


if __name__=='__main__':
    

    
    grains = [grain(-1.01,8,0.051,0),grain(1.01,8,-0.051,0),grain(-1.1,4,5,0),grain(1.1,4,-5,0),grain(-2.1,0,10,0),grain(2.1,0,-10,0),grain(-1.1,-4,50,0),grain(1.1,-4,-50,0),grain(-1.1,-8,100,0),grain(1.1,-8,-100,0)]

    T,X,Y,VX,VY,FX,FY,FpX,FpY,FvX,FvY,G = trajectoires(grains, 1)
    
    #%%
    for k in range(len(grains)//2):
        plt.plot(T,VX[2*k]**2/VX[2*k][0]**2, label = 'Vx0={}'.format(VX[2*k][0]))   
    plt.title("Ec relative en fonction du temps")
    plt.show()
    plt.legend()
    plt.grid()
    #%%
    for k in range(len(grains)//2):
        plt.scatter(VX[2*k][0],1-VX[2*k][-1]**2/VX[2*k][0]**2)   
        
    plt.title("Perte relative d'Ec en fonction de la vitesse initiale")
    

    plt.legend()
    plt.grid()
    plt.show()
    #%%
    grains = [grain(-1.01,8,0.051,0),grain(1.01,8,-0.051,0),grain(-1.1,4,5,0),grain(1.1,4,-5,0),grain(-2.1,0,10,0),grain(2.1,0,-10,0),grain(-1.1,-4,50,0),grain(1.1,-4,-50,0),grain(-1.1,-8,100,0),grain(1.1,-8,-100,0)]
    
    fig = plt.figure(1,figsize=(6,6))
    plt.axis("equal")
    ax = plt.axes(xlim=(-10,10), ylim=(-10,10))
    fenetre = fig, ax
    ani = animer(grains)
    plt.show()
    
    
    