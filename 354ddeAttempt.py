#DDE numerical solution attempt: #354
import numpy as np
import matplotlib.pyplot as plt
import math
from jitcdde import y, t, jitcdde

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
tauB=2.

#form the differential system
sys_f=[
		alphaB*math.exp(-mu*tauB)*y(1,t-tauB)-gamB*y(0),
		alphaM*(1+K1*(math.exp(-mu*tauM)*y(2,t-tauM))**n)/(K+K1*(math.exp(-mu*tauM)*y(2,tauM))**n) - gamM*y(1),
		alphaA*y(0)*L/(KL+L)-betaA*y(0)*y(2)/(KA+y(2))-gamA*y(2)
]

#initiate, insert delays
I = jitcdde(sys_f, delays=[tauB,tauM])
#arbitrary initial values for the time being
I.constant_past( [.000001,0.000001,0.001], time=0.0 )
I.integrate_blindly(tauB)

times=np.arange(1.3,200.,.2)
x=np.zeros(len(times))
i=0
while i< len(times):
	x[i]=I.integrate(i)[0]
	i=i+1
plt.plot(times,x)
plt.show()