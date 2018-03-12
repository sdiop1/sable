# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:00:52 2018

@author: LouisC

Influence de la valeur du dt sur la symetrie d'une collision
"""

# paramètres :

# grains
R = 1
m = 1

# exterieur
L=20
e=2
y_min = -500
g = 0

# paramètres des contacts
Kc = 1e4 * m/R
restitution = 0.5
Kd = restitution**2*Kc
mu = 0
Kt = 0


Kc_p = 1e6/(2*R)
restitution_p = 0
Kd_p = restitution_p**2*Kc_p

mu_p = 1
Kt_p = 0


dt = 1e-4



if __name__=='__main__':
    

    grains = [grain(-2,-0,20,50), grain(2,0,-20,50)]
    plt.title(r'Décalage du centre de gravité des deux grains, selon la valeur de $dt$')
    liste = [0.6e-2, 2e-3, 5e-4, 1e-4]
    for k in range(4):
        dt = liste[k]
        plt.subplot(2,2,k+1)
        grains = [grain(-2,-0,20,50), grain(2,0,-20,50)]
        traj = trajectoires(grains, 0.15)
        plt.plot(traj[0][0],traj[1][0], 'b', label = r'$p_{1}$')
        plt.plot(traj[0][1],traj[1][1], 'g', label = r'$p_{2}$')
        xg = (traj[0][0] + traj[0][1])/2
        yg = (traj[1][0] + traj[1][1])/2
        plt.plot(xg, yg, 'r', label = 'G')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        if k == 1:
            plt.legend(loc = 0, fontsize = 7)
        plt.grid()
        plt.title(r'$dt = ${}'.format(dt))
    plt.show()
