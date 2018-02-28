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





def perp_vec(a,v):
    '''
    Renvoie un vecteur unitaire b perpendiculaire au vecteur a, et dans le sens du vecteur v
    '''
    x,y = a
    if y != 0:
        b = np.array([1/np.sqrt(1+(x/y)**2),-x/(y*np.sqrt(1+(x/y)**2))])
    else :
        b = np.array([-y/(x*np.sqrt(1+(y/x)**2)),1/(np.sqrt(1+(y/x)**2))])
    if np.dot(b,v) > 0 : return b
    else : return -b
    




def force_grain(i,j):
    '''
    Fonction de calcul de la force due au grain i s'appliquant sur le grain j
    '''
    xi, yi = i.pos
    xj, yj = j.pos
    vxi, vyi = i.pos
    vxj, vyj = j.vit
    
    rij = np.array([xj -xi , yj-yi])
    vij = np.array([vxj - vxi, vyj - vyi])
    
    nij = rij / np.sqrt(np.dot(rij,rij))
    tij = perp_vec(nij,vij) 
    
    delta = 2*R - np.sqrt(np.dot(rij,rij)) 
    vit_tan = np.dot(vij,tij)
    
    if np.dot(nij,vij) < 0:
        return Kc*delta*nij - min((mu*Kc*delta, Kt*vit_tan))*tij
    if np.dot(nij,vij) > 0:
        return Kd*delta*nij - min((mu*Kd*delta, Kt*vit_tan))*tij          




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
    
    
    
    
    

def force_totale(i, gravite = True, paroi = True):
    '''
    Fonction renvoyant la force totale s'appliquant sur le grain i.
    gravite : détermine si le poids est pris en compte dans la force
    paroi : détermine si la force du récipient est prise en compte
    '''
    force = np.array([0.,0.],dtype='float64')
    i.detecter_voisins()
    i.detecter_paroi()
    for j in i.voisins:
        force += force_grain(j,i)
        force += force_grain(j,i)
    if paroi == True :
        force += force_paroi(i)
        force += force_paroi(i)
    if gravite == True :
        force += np.array([0.,-m*g],dtype='float64')
    return force
    

  
  
  

def mise_a_jour(remove=True):
    '''
    Fonction de mise à jour de la liste des positions/vitesses des grains par l'algorithme de Verlet à deux pas.
    1 : calculer les nouvelles vitesses/positions de tous les grains
    2 : les "nouvelles positions/vitesses" deviennent les pos/vit actuelles
        lorsque remove vaut True : si un grain sort du silo (y < y_min) on ne le prend plus en compte dans l'étude
    '''
    for i in grains:
        i.pos_n = 2 * i.pos - i.pos_a + 1/m * force_totale(i) * dt**2
        i.vit_n = (i.pos_n - i.pos_a)/(2*dt)
    for i in grains:
        i.pos, i.pos_a = i.pos_n, i.pos
        i.vit = i.vit_n
        if i.pos[1] < y_min :
            if remove == True:
                grains.remove(i)
                
                
def trajectoires(grains, t_final):
    '''
    Renvoie les trajectoires [X,Y,VX,VY,temps] qui sont tous des arrays d'un ensemble de grains initial 'grains', entre les instants 0 et t_final
    '''
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
        mise_a_jour()
        t += dt
    
    return(np.array(X),np.array(Y),np.array(VX),np.array(VY),np.array(temps))
