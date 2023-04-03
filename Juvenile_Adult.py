#Produces Figure 2

import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from numpy import log, exp, sqrt

gamma = .005
pj=.1
pa=.1
b=10
trans = 0.1
time=80
k=.8

IC = [1,1,1]

def evoLV(X,t):
    
    J = X[0]
    A = X[1]
    v = X[2]
    
    dJdt = b*(v-1)*A-pj*v**2*J-trans*J-gamma*J*(A+J)
    dAdt = trans*J-pa*A-gamma*A*(A+J)
    
    dvdt = k*(-pj*v + (b*trans - pa*pj*v + pj**2*v**3 + pj*trans*v)/(sqrt(4*b*trans*v - 4*b*trans + pa**2 - 2*pa*pj*v**2 - 2*pa*trans + pj**2*v**4 + 2*pj*trans*v**2 + trans**2)))
    
    dxvdt = np.array([dJdt, dAdt,dvdt])

    return dxvdt

intxv = np.array(IC)
pop = odeint(evoLV, intxv, range(time+1))
    
plt.figure()
plt.subplot(211)
plt.title('Juvenile-Adult Eco-Evolutionary Dynamics')
plt.plot(pop[:,0],lw=2,color='b',label='Juvenile')
plt.axhline(y=186.24228771767451,lw=1,color='b',linestyle='dashed')
plt.plot(pop[:,1],lw=2,color='r',label='Adult')
plt.axhline(y=16.7071335456545,lw=1,color='r',linestyle='dashed')
plt.ylim(0,200)
plt.xlim(0,time)
plt.legend()
plt.subplot(212)
plt.plot(pop[:,2],lw=2,color='k',label='Strategy')
plt.axhline(y=4.485014871384016,lw=1,color='k',linestyle='dashed')
plt.ylim(0,5)
plt.xlim(0,time)
plt.legend()
plt.show()