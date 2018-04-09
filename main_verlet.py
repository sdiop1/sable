#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 12:45:34 2018
@author: sdiop
Nouvelle version du main avec l'algorithme de verlet
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
import random as rd



class grain :
    '''
    Definition des méthodes et attributs des grains
    '''
    
    def __init__(self,x0,y0,vx0,vy0):
        self.pos = np.array([x0,y0]) #position à t
        self.vit = np.array([vx0,vy0]) #vitesse à t
        self.acc = np.array([0.,0.], dtype = 'float64') #accélération à t
        
        self.pos_a = np.array([x0,y0]) - np.array([vx0,vy0])*dt  # position à t-dt
        self.pos_n = np.array([x0,y0]) #position à t+dt
        self.vit_n = np.array([vx0,vy0]) #position à t+dt
        self.force = np.array([0,0]) #force_totale(i) a t 
        self.force_voisins = np.array([0.,0.])
        self.force_paroi = np.array([0.,0.])
        self.grav = np.array([0.,0.])
        
        self.voisins = [] #liste des voisins
        self.bords = [] #liste de booléens [gauche,bas,droite], True si contact, False sinon
        
    def detecter_voisins(self):
        '''
        Renvoie la liste des grains j en contact avec le grain self
        On teste tous les grains (différents de i) pour voir si ils sont en contact
        '''
        self.voisins=[]
        for k in grains:
            if ((k.pos[0]-self.pos[0])**2+(k.pos[1]-self.pos[1])**2 < 4*R**2) and k != self:
                self.voisins.append(k)
                
                
    def detecter_paroi(self):
        self.bords = []
        if self.pos[0] < -L/2 + R:
          self.bords.append(True)
        else:
           self.bords.append(False)
        if self.pos[1] < R and (self.pos[0] < -e or self.pos[0] > e):
           self.bords.append(True)
        else:
           self.bords.append(False)
        if self.pos[0] > L/2 - R:
           self.bords.append(True)
        else:
           self.bords.append(False)



def force_grain(i,j):
    '''
    Fonction de calcul de la force due au grain i s'appliquant sur le grain j
    '''
    xi, yi = i.pos
    xj, yj = j.pos
    vxi, vyi = i.pos
    vxj, vyj = j.vit
    
    rij = np.array([xj - xi, yj - yi])
    x, y = rij
    vij = np.array([vxj - vxi, vyj - vyi])
    r = np.linalg.norm(rij)
    
    u_r = 1/r * np.array([x, y])
    u_t = 1/r * np.array([-y, x])
    
    v_t = np.dot(vij, u_t)
    v_r = np.dot(vij, u_r)
    
    if np.dot(vij,u_r) < 0 :
        K = Kc
    else :
        K = Kd
    
    F_r = K*(2*R - r)
    if np.abs(Kt * v_t) < np.abs(mu * F_r) : 
        F_t = - Kt * v_t
    else :
        F_t = - mu * np.abs(F_r) * np.sign(v_t)
        
    return F_r * u_r + F_t * u_t
       




def force_paroi(i):
    '''
    Fonction de calcul de la force s'appliquant sur le grain i à cause de la paroi
    '''    
    x,y = i.pos
    vx,vy = i.vit
    f = np.array([0.,0.])
    
    if i.bords[0] == True :
        alpha = -(x-R+L/2)
        if vx > 0 :
            f += np.array([Kd_p*alpha,-min((mu_p*Kc_p*alpha, Kt_p*vy))])
        if vx <= 0 :
            f += np.array([Kc_p*alpha,-min((mu_p*Kd_p*alpha, Kt_p*vy))])
    if i.bords[2] == True :
        alpha = x+R-L/2
        if vx > 0 :
            f += np.array([-Kc_p*alpha,-min((mu_p*Kc_p*alpha, Kt_p*vy))])
        if vx <= 0 :
            f += np.array([-Kd_p*alpha,-min((mu_p*Kd_p*alpha, Kt_p*vy))])
    if i.bords[1] == True :
        alpha = R-y
        if vy > 0 :
            f += np.array([-min((mu_p*Kc_p*alpha, Kt_p*vx)), Kd_p*alpha])
        if vy <= 0 :
            f += np.array([-min((mu_p*Kc_p*alpha, Kt_p*vx)), Kc_p*alpha])
    return f
    
    
    
    
    

def force_totale(i):
    '''
    Fonction renvoyant la force totale s'appliquant sur le grain i.
    gravite : détermine si le poids est pris en compte dans la force
    paroi : détermine si la force du récipient est prise en compte
    '''
    force = np.array([0.,0.],dtype='float64')
    i.detecter_voisins()
    i.detecter_paroi()
    
    for j in i.voisins:
        fg = force_grain(j,i) 
        force += fg
        
    if paroi == True :
        fp = force_paroi(i)  
        force += fp 
        
    if gravite == True :
        fgrav = np.array([0.,-m*g],dtype='float64')
        force += fgrav
        
    return force
    

  
  
  

def mise_a_jour(remove=True):
    '''
    Fonction de mise à jour de la liste des positions/vitesses des grains par l'algorithme de Verlet à deux pas.
    1 : calculer les nouvelles vitesses/positions de tous les grains
    2 : les "nouvelles positions/vitesses" deviennent les pos/vit actuelles
        lorsque remove vaut True : si un grain sort du silo (y < y_min) on ne le prend plus en compte dans l'étude
    '''
    for i in grains:
        f = force_totale(i)
        i.force = f 
        i.pos_n = 2 * i.pos - i.pos_a + 1/m * f * dt**2
        i.vit_n = (i.pos_n - i.pos_a)/(2*dt)
    for i in grains:
        i.pos, i.pos_a = i.pos_n, i.pos
        i.vit = i.vit_n
        if remove == True:
            if i.pos[1] < y_min :
                grains.remove(i)
                
                
                
def trajectoires(grains, t_final):
    '''
    Renvoie les trajectoires [T,X,Y,VX,VY] qui sont tous des arrays d'un ensemble de grains initial 'grains', entre les instants 0 et t_final
    '''
    t = 0
    X = [[] for k in range(len(grains))]
    Y = [[] for k in range(len(grains))]
    VX = [[] for k in range(len(grains))]
    VY = [[] for k in range(len(grains))]

    
    temps = []

    while t < t_final:
        for i in range(len(grains)) :
            X[i].append(grains[i].pos[0])
            Y[i].append(grains[i].pos[1])
            VX[i].append(grains[i].vit[0])
            VY[i].append(grains[i].vit[1])
        temps.append(t)
        mise_a_jour()
        t += dt
        
    for k in range(len(X)): X[k] = np.array(X[k])
    for k in range(len(Y)): Y[k] = np.array(Y[k])
    for k in range(len(VX)): VX[k] = np.array(VX[k])
    for k in range(len(VY)): VY[k] = np.array(VY[k])
    return(temps, X, Y, VX, VY)