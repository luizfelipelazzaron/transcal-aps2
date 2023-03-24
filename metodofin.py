import numpy as np
import matplotlib.pyplot as plt
import math

# matriz3d = np.zeros((6,100))
# matriz3d[0,:,:] = 100

# alpha = 0.25
# dx = 0.1
# dt = 0.01
# for k in np.arange(0,99):
#     for i in np.arange(1,5):
#         for j in np.arange(1,5):
#             matriz3d[i,j,k+1] = matriz3d[i,j,k] + (alpha*dt/dx**2)*(matriz3d[i,j+1,k] - 2*matriz3d[i,j,k] + matriz3d[i,j-1,k] + matriz3d[i+1,j,k] - 2*matriz3d[i,j,k] + matriz3d[i-1,j,k])
# # print(matriz3d)
# print(matriz3d[1,1,99])
# plt.imshow(matriz3d[:,:,99], cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()

def metodo_finito(rho, cp, k, h, Tinf, Tb, raio, l, dx, Ttot, dt):
    """
    Método de diferenças finitas para resolução de problemas de transferência de calor 1D
    """
    # Definindo o número de nós
    N = int(l/dx)
    # Definindo o número de pontos de tempo
    M = int(Ttot/dt)
    # Definindo o coeficiente de difusão
    alpha = k/(rho*cp)
    # Definindo o vetor de posições
    x = np.arange(0,l+dx+1,dx)
    # Definindo o vetor de tempos
    t = np.arange(0,Ttot+dt,dt)
    # Definindo a matriz de temperaturas
    T = np.zeros((N,M))
    # Definindo a condição de contorno
    T[:,:] = Tinf
    T[0,:] = Tb
    T[N-1,:] = 25 + 273

    for k in np.arange(0,M-1):
        for i in np.arange(1,N-1):
            T[i,k+1] = T[i,k] + (alpha*dt)*((1/dx**2)*(T[i+1,k] - 2*T[i,k] + T[i-1,k]) + ((h*2*dx/(raio*k))*(Tinf - T[i,k])))
    return T

# Definindo os parâmetros

rho = 2700
cp = 896
k = 180
h = 50
Tinf = 50 + 273
Tb = 100 + 273
raio = 5e-3/2
l = 300e-3
dx = 0.001
P = 2*math.pi*raio
A = math.pi*raio**2
alpha = k/(rho*cp)
dt = 0.9*dx**2/(alpha*((h*P/(A*k))+2))
Ttot = 500
Tl = 25 + 273

Temperaturas = metodo_finito(rho, cp, k, h, Tinf, Tb, raio, l, dx, Ttot, dt)

def temp_final(Tb, Tinf, Tl, P, A, h, K, L):
    Tf = []
    # Diferença de temperatura entre a base e a temperatura infinita (Kelvin)
    teta_b = Tb - Tinf
    teta_l = Tl - Tinf
    #Parâmetros
    M = teta_b*(h * P * K * A)**0.5
    m = ((h * P) / (K * A))**0.5
    # Posições ao longo do comprimento (m)
    tamanhos = np.linspace(0, L, 1000)
    for x in tamanhos:
      #Diferença de temperatura
      teta = teta_b * (((teta_l/teta_b)*np.sinh(m*x) + np.sinh(m*(L-x))) / (np.sinh(m*L)))
      
      # Temperaturas em cada posição ao longo do comprimento (Kelvin)
      Tf.append(Tinf + teta)
    return Tf, tamanhos

# Tf, tamanhos = temp_final(Tb, Tinf, Tl, P, A, h, k, l)
print(Temperaturas)
# Plotando o gráfico
# plt.plot(tamanhos, Tf, label='Solução analítica')
# plt.xlabel('Posição (m)')
# plt.ylabel('Temperatura (K)')
# plt.title('Temperatura ao longo do comprimento')
# plt.grid(True)
# plt.show()
