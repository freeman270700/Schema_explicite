#Réponse à la question
# On voit que la convergence ce fait à partir de 100 itérations


import matplotlib.pyplot as plt
import numpy as np
from numpy import nbytes, zeros,array,dot,linspace, linalg
from math import *
import matplotlib.pyplot as plt

start_time = time.time()
K = 2*10**-6
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
# Graphe de Convergence

# Graphe de Convergence
plt.figure(1)
plt.plot(np.log(err))
plt.title("Graphe de Convergence")

# Variation de temperature (imagesc)
plt.figure(2)
plt.imshow(u, extent=[0, 12, 0, 100], aspect='auto')
plt.title("Variation de temperature (imagesc)")
plt.colorbar()

# Variation de Temperature(contourf)
plt.figure(3)
profondeur = np.arange(0, Nz*dz+dz, dz)
plt.contourf(temps, -profondeur, u)
plt.title("Variation de Temperature(contourf)")
plt.colorbar()

# Tracé de la température en fonction du temps à différentes profondeurs
plt.figure(4)
plt.plot(temps, u[0,:], label="0m")
plt.plot(temps, u[20,:], label="5m")
plt.plot(temps, u[40,:], label="10m")
plt.plot(temps, u[60,:], label="15m")
plt.plot(temps, u[80,:], label="20m")
plt.xlabel("Temps (mois)")
plt.ylabel("Temperature (C)")
plt.legend()

plt.show()


end_time = time.time()
# Calcul du temps d'exécution en secondes
elapsed_time = end_time - start_time
print(f"Temps d'exécution : {elapsed_time:.2f} secondes")


