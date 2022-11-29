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
temp = np.ones((Nz+1, Nt+1))
pas = 12/Nt
u = 15 * temp
for i in range(0, Nt+1):
    u[0, i] = 15 - 10 * sin(2*pi*i*pas/12)
maxiter = 500
for iter in range(0 , maxiter) : 
    uold = u
    u[:,0] = uold[:,-1]
    for i in range(1, Nt+1):
        profondeur = (u[0:(Nz)-2,i-1]) 
        profondeur=- 2*(u[1:(Nz)-1,i-1])
        profondeur=+ (u[2:(Nz),i-1])
        temps_1D = K*profondeur  
        u[1:Nz,i] =+ u[1:Nz,i-1] 
        u[2:Nz,i]=+ temps_1D
        u[-1,i] = u[len(u)-1,i]
    #trouver le maximum en valeur absolue entre deux solutions
    diff = abs(uold - u)
    diffmax = diff.max()
    print(diffmax)
    if diffmax < 1e-4 : 
        break
    end
#afficher valeur err
print(diffmax)
#sortie graphique
