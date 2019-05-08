#Solve for the fixed points with varied L concentration
import numpy as np
import matplotlib.pyplot as plt
import math

#define parameters
mu = 3.03*(10**(-2))
n = 2.0
tauM = 0.1
K = 7200
K1 = .0252

betaA = 2.15*(10**4)
alphaB = .0166
alphaM = .997
alphaA = 1.76*(10**4)
KL=970
KA=1950
L=50
L2=45
L3=40
gamM1 = 0.411
gamA1 = 1.35*(10**(-2))
gamB1 = 8.33*(10**(-4))
gamM = gamM1+mu
gamA = gamA1+mu
gamB = gamB1+mu
tauB=2

#initiate arrays to hold values
A = np.arange(0., 40., 1.)
LHS=np.zeros(len(A))
RHS=np.zeros(len(A))
RHS2=np.zeros(len(A))
RHS3=np.zeros(len(A))

#Solve LHS, RHS(s)
i=0
while i < len(A):
    LHS[i] = (1+K1*(math.exp(-mu*tauM)*A[i])**n)/(K+K1*(math.exp(-mu*tauM)*A[i])**n)

    h = L / (KL + L)
    h2 = L2 / (KL + L2)
    h3 = L3 / (KL + L3)
    g = A[i] / (KA + A[i])
    X = (gamM * gamB * gamA * math.exp(mu * tauB)) / (alphaM * alphaB * alphaA)
    RHS[i] = X * A[i] / (h - (betaA * g / alphaA))
    RHS2[i] = X * A[i] / (h2 - (betaA * g / alphaA))
    RHS3[i] = X * A[i] / (h3 - (betaA * g / alphaA))
    i=i+1

#Ploy
plt.plot(A,LHS,linestyle='--',label='LHS')
plt.plot(A,RHS,label='L=50 (\u03bcM)')
plt.plot(A,RHS2,label='L=45 (\u03bcM)')
plt.plot(A,RHS3,label='L=40 (\u03bcM)')
plt.title("Steady State Values")
plt.xlabel('A Concentration (\u03bcM)')
plt.legend(loc='upper left')
plt.show()


