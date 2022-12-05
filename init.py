import matplotlib.pyplot as plt
import numpy as np
from numpy import nbytes, zeros,array,dot,linspace, linalg
from math import *
import matplotlib.pyplot as plt

K = 2*10e-6
dz = 0.25
Nz = 400
Nt = 5000
dt = (365*24*60*60)/Nt
u = np.ones((Nz+1, Nt+1))
#affiche la taille de temp
temp = np.linspace(0,12,Nt+1)
u = 15 * u
#pour i allant de 0 à Nt
for i in range(0,Nt+1):
    u[0, i] = 15 - 10 * sin(2*pi*temp[i]/12)
maxiter = 500
for iter in range(0 , maxiter) : 
    uold = u.copy()
    u[:,0] = uold[:,Nt]
    for i in range(1, Nt+1):
        #profondeur[i] = ((u[0:(Nz-2),i-1]) - 2*(u[1:(Nz-1),i-1]) + (u[2:(Nz),i-1]))/dz**2
        temps_1D = K*profondeur  
        u[1:Nz-2,i] = dt*temps_1D[1:Nz-2] + u[1:Nz-2,i-1] 
        u[Nz-1,i] = u[Nz-2,i]
        
    #trouver le maximum en valeur absolue entre deux solutions u et uold
    diff = np.max(np.abs(u-uold))
    if(diff < 1e-4):
        break

#afficher valeur err

#sortie graphique variation de la température en fonction de la profondeur
plt.figure()
#mettre les abscisses en fcontion du temps et les ordonnées en fonction de la temperature

print(u[0:(Nz-2),0])
#mettre les ordonnées entre 0 et 15
plt.ylabel('Température (°C)')
plt.plot(temp,u[0,:], 'r')
plt.plot(temp,u[20,:], 'r')
plt.plot(temp,u[40,:], 'b')
plt.plot(temp,u[60,:], 'b')
plt.plot(temp,u[80,:], 'g')
plt.show()


