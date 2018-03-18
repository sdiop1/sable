#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 12:45:34 2018
@author: sdiop

Code principal avec :
    - algorithme de Verlet
    - modèle de dissipation
    - degré de liberté de rotation
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
        self.theta = 0.
        self.omega = 0
        
        self.theta_n = 0.
        self.omega_n = 0.
        self.theta_a = 0
        self.pos_a = np.array([x0,y0]) - np.array([vx0,vy0])*dt  # position à t-dt
        self.pos_n = np.array([x0,y0]) #position à t+dt
        self.vit_n = np.array([vx0,vy0]) #position à t+dt
        
        self.force = np.array([0.,0.]) #force_totale(i) a t 
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





def perp_vec(a,v):
    '''
    Renvoie un vecteur unitaire b perpendiculaire au vecteur a, et dans le sens du vecteur v
    '''
    X=a[0]
    Y=a[1]
    if Y != 0:
        b = np.array([1/np.sqrt(1+(X/Y)**2),-X/(Y*np.sqrt(1+(X/Y)**2))])
    else :
        b = np.array([-Y/(X*np.sqrt(1+(Y/X)**2)),1/(np.sqrt(1+(Y/X)**2))])
    if scalaire(b,v) > 0 : return b
    else : return -b
    
    
def scalaire(u,v):
    n_u = np.shape(u)[0]
    n_v = np.shape(v)[0]
    s = 0
    for k in range(min(n_u,n_v)):
        s += u[k] * v[k]
    return s

def gv(a,seuil):
    return(1-np.exp(-a/seuil))


def force_grain(i,j):
    '''
    Fonction de calcul de la force due au grain i s'appliquant sur le grain j
    '''
    xi, yi = i.pos
    xj, yj = j.pos
    vxi, vyi = i.pos
    vxj, vyj = j.vit
    omi = i.omega
    omj = j.omega
    
    rij = np.array([xj -xi , yj-yi])
    rij3d = np.array([xj -xi , yj-yi, 0])
    vij = np.array([vxj - vxi, vyj - vyi])
    vij3d = np.array([vxj - vxi, vyj - vyi, 0])
    omij = omj - omi
    
    nij = rij / np.sqrt(np.dot(rij,rij))
    nij3d = rij3d / np.sqrt(np.dot(rij,rij))
    tij = perp_vec(nij,vij)
    
    delta = 2*R - np.sqrt(np.dot(rij,rij)) 
    vit_tan = scalaire(vij3d + omij*(R-delta/2)*np.cross(nij3d, np.array([0.,0.,1.])) , tij)
    vit_norm = np.dot(vij,nij)
    
    force_normale = Kc*delta - gamma*gv(delta,zeta)*vit_norm
    force_tangent = - min((mu*np.abs(force_normale), Kt*np.abs(vit_tan))) 
    moment = -np.cross( (R-delta/2)*nij , force_tangent*tij )

    return force_normale*nij + force_tangent*tij , moment
    




def force_paroi(i):
    '''
    Fonction de calcul de la force s'appliquant sur le grain i à cause de la paroi
    '''    
    x,y = i.pos
    vx,vy = i.vit
    om = 0
    
    force = np.array([0.,0.])
    moment = 0.
    
    
    if i.bords[0] == True :
        delta = -(x-R+L/2)
        force_normale = Kc_p*delta - gamma_p*gv(delta,zeta_p)*vx
        force_tangent = -min((mu_p*np.abs(force_normale), Kt_p*np.abs(vy+(R-delta/2)*om)))*np.sign(vy)
        force += np.array([force_normale,force_tangent])
        moment += -np.cross( (R-delta/2)*np.array([-1.,0.]) , force_tangent*np.array([0.,1.]) )
        
    if i.bords[2] == True :
        delta = x+R-L/2
        force_normale = - Kc_p*delta - gamma_p*gv(delta,zeta_p)*vx
        force_tangent = -min((mu_p*np.abs(force_normale), Kt_p*np.abs(vy-(R-delta/2)*om)))*np.sign(vy)
        force += np.array([force_normale,force_tangent])
        moment += -np.cross( (R-delta/2)*np.array([1.,0.]) , force_tangent*np.array([0.,1.]) )

    if i.bords[1] == True :
        delta = R-y
        force_normale = Kc_p*delta - gamma_p*gv(delta,zeta_p)*vy
        force_tangent = -min((mu_p*np.abs(force_normale), Kt_p*np.abs(vx+(R-delta/2)*om)))*np.sign(vx)
        force += np.array([force_tangent, force_normale])
        moment += -np.cross( (R-delta/2)*np.array([0.,-1.]) , force_tangent*np.array([1.,0.]) )
        
    return force, -moment
    
    
    
    
    

def force_totale(i):
    '''
    Fonction renvoyant la force totale s'appliquant sur le grain i.
    gravite : détermine si le poids est pris en compte dans la force
    paroi : détermine si la force du récipient est prise en compte
    '''
    force = np.array([0.,0.],dtype='float64')
    moment = 0
    i.detecter_voisins()
    i.detecter_paroi()
    i.force_voisins = np.array([0.,0.])
    i.force_paroi = np.array([0.,0.])
    i.grav = np.array([0.,0.])
    for j in i.voisins:
        fg = force_grain(j,i)[0]
        moment += force_grain(j,i)[1]
        force += fg
        i.force_voisins += fg 
    if paroi == True :
        fp = force_paroi(i)[0]
        moment += force_paroi(i)[1]
        force += fp 
        i.force_paroi += fp
    if gravite == True :
        fgrav = np.array([0.,-m*g],dtype='float64')
        force += fgrav
        i.grav = fgrav
    i.force += force
    
    return force, moment
    

  
  
  

def mise_a_jour():
    '''
    Fonction de mise à jour de la liste des positions/vitesses des grains par l'algorithme de Verlet à deux pas.
    1 : calculer les nouvelles vitesses/positions de tous les grains
    2 : les "nouvelles positions/vitesses" deviennent les pos/vit actuelles
        lorsque remove vaut True : si un grain sort du silo (y < y_min) on ne le prend plus en compte dans l'étude
    '''
    for i in grains:
        f,mom = force_totale(i)
        i.force = f 
        i.pos_n = 2 * i.pos - i.pos_a + 1/m * f * dt**2
        i.vit_n = (i.pos_n - i.pos_a)/(2*dt)
        i.theta_n = 2 * i.theta - i.theta_a + mom * dt**2
        i.omega_n = (i.theta_n - i.theta_a)/(2*dt)
    for i in grains:
        i.pos, i.pos_a = i.pos_n, i.pos
        i.theta, i.theta_a = i.theta_n *0, i.theta *0
        i.vit = i.vit_n
        i.omega = i.omega_n *0
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
    FX = [[] for k in range(len(grains))] #forces totales en x sur les grains
    FY = [[] for k in range(len(grains))]  
    FpX = [[] for k in range(len(grains))]#forces de la paroi en y sur les grains 
    FpY = [[] for k in range(len(grains))]
    FvX = [[] for k in range(len(grains))]#forces des voisisns en x sur les grains
    FvY = [[] for k in range(len(grains))]
    G = [[] for k in range(len(grains))]
    
    temps = [0]
    for i in range(len(grains)) :
        X[i].append(grains[i].pos[0])
        Y[i].append(grains[i].pos[1])
        VX[i].append(grains[i].vit[0])
        VY[i].append(grains[i].vit[1])
        FX[i].append(grains[i].force[0])
        FY[i].append(grains[i].force[1])
        FpX[i].append(grains[i].force_paroi[0])
        FpY[i].append(grains[i].force_paroi[1])
        FvX[i].append(grains[i].force_voisins[0])
        FvY[i].append(grains[i].force_voisins[1])
        G[i].append(grains[i].grav[1])
    
    while t < t_final:
        for i in range(len(grains)) :
            X[i].append(grains[i].pos[0])
            Y[i].append(grains[i].pos[1])
            VX[i].append(grains[i].vit[0])
            VY[i].append(grains[i].vit[1])
            FX[i].append(grains[i].force[0])
            FY[i].append(grains[i].force[1])
            FpX[i].append(grains[i].force_paroi[0])
            FpY[i].append(grains[i].force_paroi[1])
            FvX[i].append(grains[i].force_voisins[0])
            FvY[i].append(grains[i].force_voisins[1])
            G[i].append(grains[i].grav[1])
        temps.append(t)
        mise_a_jour()
        t += dt
    return(np.array(temps),np.array(X),np.array(Y),np.array(VX),np.array(VY),np.array(FX),np.array(FY),np.array(FpX),np.array(FpY),np.array(FvX),np.array(FvY),np.array(G))
