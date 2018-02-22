# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:00:52 2018

@author: LouisC
"""
import matplotlib.pyplot as plt 


# paramètres :

# grains
R = 1
m = 1

# exterieur
L=200
e=0
y_min = -500
g = 10

# paramètres des contacts
Kc = 1e4 * m/R
restitution = 0.7
Kd = restitution**2*Kc
mu = 0.5
#Kt = 10


Kc_p = 1e6/(2*R)
restitution_p = 0.7
Kd_p = restitution_p**2*Kc_p

mu_p = 0.5
#Kt_p = Kt


dt = 1e-4



if __name__=='__main__':
    

    for k in range(4):

        Kt = 10*(k+1)
        Kt_p = Kt

        grains = [grain(0,1.,1,0)]
        #plt.title(r'Evolution du grain pour Kt = Kt_p =.{}'.format(Kt))
        
        traj = trajectoires(grains, 1)
        plt.subplot(2,2,1)
        plt.plot(traj[4],traj[2][0], label = r'$Kt=${}'.format(Kt))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$Vx$')
        plt.title(r'$V_{x}(t)$')
        #plt.ylim(max(traj[2][0]) + max(traj[2][0])/2, ymin = min(traj[2][0]) - min(traj[2][0])/2 )
        
        
        plt.subplot(2,2,2)
        plt.plot(traj[4],traj[0][0], label = r'$Kt=${}'.format(Kt))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$x$')
        plt.title(r'$x(t)$')
    
        plt.subplot(2,2,3)
        plt.plot(traj[0][0],traj[1][0], label = r'$Kt=${}'.format(Kt))
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        plt.title(r'$y(x)$')
        
        plt.legend(loc = 0, fontsize = 7)
        plt.grid()
        
        plt.show()
        
######CCL : A gde vitesse pas de changement, a faible vitesse ---> Kt = 10 est plus physique 