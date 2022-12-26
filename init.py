#Réponse à la question
# On voit que la convergence ce fait à partir de 100 itérations


import matplotlib.pyplot as plt
import numpy as np
from numpy import nbytes, zeros,array,dot,linspace, linalg
from math import *
import matplotlib.pyplot as plt
import time

start_time = time.time()
K = 2*10**-6
dz = 0.25
Nz = 400
Nt = 5000
dt = (365*24*3600)/Nt
temps = np.linspace(0, 12, Nt+1)
u = 15 * np.ones((Nz+1, Nt+1))
u[0, :] = 15 - 10*np.sin(2*np.pi*temps/12)
maxiter = 500
err = np.ones(maxiter) # Initialise les valeurs de erreurs à 1

for iter in range(maxiter):
    uold = u.copy()  # Sauvegarde la valeur de u dans uold
    u[:, 0] = uold[:, -1]
    for i in range(1, Nt+1):
        profondeur = (u[0:len(u)-2, i-1]-2*u[1:len(u)-1, i-1]+u[2:len(u), i-1])/dz**2
        temps_1D = K*profondeur
        u[1:len(u)-1, i] = temps_1D*dt + u[1:len(u)-1, i-1]
        u[-1, i] = u[-2, i]
        
    err[iter] = np.max(np.abs(u-uold))  # Calcul l'erreur absolue entre u et uold
    if err[iter] < 1E-4:
        break


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
print(f"Temps d'execution : {elapsed_time:.2f} secondes")


