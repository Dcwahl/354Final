import numpy as np
import matplotlib.pyplot as plt
import math

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
A = np.arange(2.3, 70., .1)

g = A / (KA + A)
X = (gamM * gamB * gamA * math.exp(mu * tauB)) / (alphaM * alphaB * alphaA)
f =(1+K1*(math.exp(-mu*tauM)*A)**n)/(K+K1*(math.exp(-mu*tauM)*A)**n)
m=(1/f)*(X*A+f*(betaA/alphaA)*g)
L=(KL*m)/(1-m)

print(L)
i=0
while i <len(A):
    A[i]=math.log(A[i])
    i=i+1
plt.plot(L,A)
plt.title("3 Variable Steady States")
plt.xlabel('L Concentration (\u03bcM)')
plt.ylabel('A Steady State Concentration (\u03bcM)')
plt.axvline(x=39, color='k', linestyle='--')
plt.axvline(x=55.5, color='k', linestyle='--')
plt.show()