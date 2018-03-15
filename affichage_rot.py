#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
import matplotlib.transforms as trf 
from random import random

# paramètres :

dt = 1e-4





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
    if frame_nb%1 == 0:
        for artist in plt.gca().lines + plt.gca().collections:
            artist.remove()
        patches = []
        for i in grains :
            x = i.pos[0]
            y = i.pos[1]
            circle = pat.Circle((x, y), R)
            trans = trf.Affine2D().rotate_around(x,y,i.thet%(2*np.pi)) 
            circle.set_transform(trans)
            patches.append(circle)
        p = PatchCollection(patches, alpha=0.5, linestyle = 'dotted',linewidth = 3, color = 'r')
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

