#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 17:38:47 2018

@author: sdiop
"""

import main_verlet
import main_euler

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


dt = 1e-5



paroi = False




def mise_a_jour_euler():
    # 1 : calculer les nouvelles vitesses/positions de tous les grains
    # 2 : les "nouvelles positions/vitesses" deviennent les pos/vit actuelles
    #     si un grain sort du silo on ne le prend plus en compte dans l'étude
    for i in grains:
        i.acc = 1/m * force_totale(i, False)
        i.vit_n = i.vit + dt * i.acc
        i.pos_n = i.pos + dt * i.vit_n
    for i in grains:
        i.pos = i.pos_n
        i.vit = i.vit_n
        if i.pos[1] < y_min :
            grains.remove(i)
            
            
            
def mise_a_jour_verlet(remove=True):
    global R
    # 1 : calculer les nouvelles vitesses/positions de tous les grains, algo de Verlet
    # 2 : les "nouvelles positions/vitesses" deviennent les pos/vit actuelles
    #     si un grain sort du silo on ne le prend plus en compte dans l'étude
    for i in grains:
        # i.vit_n = i.vit + dt * i.acc
        # i.pos_n = i.pos + dt * i.vit_n
        i.acc = 1/m * force_totale(i, False)
        i.pos_n = i.pos + i.vit * dt + i.acc * dt**2 *0.5
        i.acc = 1/m * force_totale(i, True)
        i.vit_n = i.vit + i.acc * dt/2
        i.vit_n = i.vit_n + i.acc * dt/2
    for i in grains:
        i.pos = i.pos_n
        i.vit = i.vit_n
        if i.pos[1] < y_min :
            if remove == True:
                grains.remove(i)
                
                
def trajectoires_euler(grains, t_final):
    #renvoie les trajectoires [X,Y,VX,VY,temps] qui sont tous des arrays d'un ensemble de grains initial 'grains', entre les instants 0 et t_final
    t = 0
    X = [[] for k in range(len(grains))]
    Y = [[] for k in range(len(grains))]
    VX = [[] for k in range(len(grains))]
    VY = [[] for k in range(len(grains))]
    temps = []
    while t < t_final:
        temps.append(t)
        for i in range(len(grains)) :
            X[i].append(grains[i].pos[0])
            Y[i].append(grains[i].pos[1])
            VX[i].append(grains[i].vit[0])
            VY[i].append(grains[i].vit[1])
        mise_a_jour_euler()
        t += dt
    
    return(np.array(X),np.array(Y),np.array(VX),np.array(VY),np.array(temps))
    

def trajectoires_verlet(grains, t_final):
    #renvoie les trajectoires [X,Y,VX,VY,temps] qui sont tous des arrays d'un ensemble de grains initial 'grains', entre les instants 0 et t_final
    t = 0
    X = [[] for k in range(len(grains))]
    Y = [[] for k in range(len(grains))]
    VX = [[] for k in range(len(grains))]
    VY = [[] for k in range(len(grains))]
    temps = []
    while t < t_final:
        temps.append(t)
        for i in range(len(grains)) :
            X[i].append(grains[i].pos[0])
            Y[i].append(grains[i].pos[1])
            VX[i].append(grains[i].vit[0])
            VY[i].append(grains[i].vit[1])
        mise_a_jour_verlet()
        t += dt
    
    return(np.array(X),np.array(Y),np.array(VX),np.array(VY),np.array(temps))
    
    
    
    
    
    
if __name__=='__main__':
    
    grains = [grain(-2,-0,20,0), grain(2,0,-20,0)]
    Xv,Yv,VXv,VYv,T = trajectoires_verlet(grains,1)
    
    grains = [grain(-2,-0,20,0), grain(2,0,-20,0)]
    Xe,Ye,VXe,VYe,T  = trajectoires_euler(grains,1)
    
    plt.plot(T, (Xe[0]+Xe[1])/2, 'r', label = 'euler')
    plt.plot(T, (Xv[0]+Xv[1])/2, 'b', label = 'verlet')
    
    plt.legend()
    
    plt.show()
    
    
    
    
    