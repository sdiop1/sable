# -*- coding: utf-8 -*-
"""
Code avec la methode de verlet : main 
"""

#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
import random as rd

#imports de 'C:\\Users\\diops\\Documents\\Programmation\\sable'

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





class grain :
    
    def __init__(self,x0,y0,vx0,vy0):
        self.pos = np.array([x0,y0]) #position initiale du centre
        self.vit = np.array([vx0,vy0]) #vitesse initiale nulle
        self.acc = np.array([0,0],dtype='float64') #accélération
        self.pos_n = np.array([x0,y0])
        self.vit_n = np.array([vx0,vy0]) 
        self.color = (rd.randint(0,1),rd.randint(0,1),rd.randint(0,1),rd.randint(0,1))
         
        self.voisins = [] #liste des voisins
        self.bords = [] #liste de booléens [gauche,bas,droite], True si contact, False sinon
        
    def detecter_voisins(self):
        # on teste tous les grains (différents de i) pour voir si ils sont en contact
        self.voisins=[]
        for k in grains:
            if ((k.pos[0]-self.pos[0])**2+(k.pos[1]-self.pos[1])**2 < 4*R**2) and k != self:
                self.voisins.append(k)
                
                
    def contact_paroi(self):
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
    x,y = a
    if y != 0:
        b = np.array([1/np.sqrt(1+(x/y)**2),-x/(y*np.sqrt(1+(x/y)**2))])
    else :
        b = np.array([-y/(x*np.sqrt(1+(y/x)**2)),1/(np.sqrt(1+(y/x)**2))])
    if np.dot(b,v) > 0 : return b
    else : return -b
    




def force_grain(i,j, nouveau = False):
    '''A partir de deux grains i et j en contact renvoie FTo = np.array([FTox, FToy])  i-->j'''
    if nouveau == True :
        xi, yi = i.pos_n
        xj, yj = j.pos_n
        vxi, vyi = i.pos_n
        vxj, vyj = j.vit_n
    if nouveau == False :
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




def force_paroi(i, nouveau = False):
    if nouveau == True :
        x,y = i.pos_n
        vx,vy = i.vit_n
    if nouveau == False :
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
    
    
    
    
    

def force_totale(i, nouveau = False):
    force = np.array([0,-m*g],dtype='float64')
    i.detecter_voisins()
    i.contact_paroi()
    for j in i.voisins:
        if nouveau == True:
            force += force_grain(j,i, True)
        else :
            force += force_grain(j,i, False)
    if paroi == True :
        if nouveau == True :
            force += force_paroi(i, True)
        else :
            force += force_paroi(i, False)
    return force
    

  
  
  

def mise_a_jour(remove=True):
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




def afficher(grains0):
    global fenetre
    fig, ax = fenetre
    patches = []
    for i in grains0 :
        x = i.pos[0]
        y = i.pos[1]
        circle = pat.Circle((x, y), R)
        patches.append(circle)
    p = PatchCollection(patches, alpha=0.5, linestyle = 'solid', color = 'k')
    ax.add_collection(p)
    plt.show()


def update(frame_nb):
    global fenetre
    for artist in plt.gca().lines + plt.gca().collections:
        artist.remove()
    patches = []
    for i in grains :
        x = i.pos[0]
        y = i.pos[1]
        circle = pat.Circle((x, y), R)
        patches.append(circle)
    p = PatchCollection(patches, alpha=0.5, linestyle = 'solid', color = 'r')
    ax.add_collection(p)
    plt.title('t = {:2f}'.format(frame_nb*dt)) 
    mise_a_jour()
    



def animer(grains0, film = False):
    global fenetre
    #représente l'évolution temporelle et spatiale d'un ensemble initial de grains (grains0)
    grains = grains0
    ani = animation.FuncAnimation(fig, update,  blit=False, interval=1, repeat=True)
    return ani
  
def trajectoires(grains, t_final):
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
        mise_a_jour()
        t += dt
    
    return(np.array(X),np.array(Y),np.array(VX),np.array(VY),np.array(temps))