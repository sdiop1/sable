#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
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

  
