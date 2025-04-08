import control as ct 
import numpy as np 
import matplotlib.pyplot as plt 

s = ct.tf('s')
sol = 4.3/(2* s *(1+0.5* s )*(1+0.25* s ))
zeta = 0.2
omega_n = 12
s1 = -zeta * omega_n + omega_n * np.sqrt (1- zeta**2)*1j
g = 4.3/(2*s1*(1+0.5*s1)*(1+0.25*s1))
varphi = np.angle (s1)
varphi_c = np.pi - np.angle ( g )
theta_p =( varphi - varphi_c )/2.0
theta_c =( varphi + varphi_c )/2.0
pc = -zeta * omega_n - omega_n * np.sqrt (1- zeta **2)/ np.tan ( theta_p )
s1 = -zeta * omega_n + omega_n * np.sqrt (1- zeta **2)*1j
g =4.3/(2*s1*(1+0.5*s1)*(1+0.25*s1))
varphi = np.angle (s1)
varphi_c = np.pi - np.angle(g)
theta_p =(varphi - varphi_c)/2.0
theta_c =(varphi + varphi_c)/2.0
pc = -zeta * omega_n - omega_n * np.sqrt (1- zeta **2)/ np.tan ( theta_p )
zc = -zeta * omega_n - omega_n * np.sqrt (1- zeta **2)/ np.tan ( theta_c )
print(pc, zc)
alpha = zc / pc 
kc = np.absolute (s1-pc ) / (np.absolute(g)* np.absolute(s1-zc)) # magn :

lead_compensator = ct.zpk([zc], [pc], kc)

scl = ct.feedback ( sol )# without compensation 
sys = ct.feedback ( sol * lead_compensator )

fig, ax = plt.subplots ()
t = np.arange (0, 10.01, 0.01)
T, y = ct.step_response (scl, t )
ax.plot ( T, y, linewidth=4, label ='w.o. compensation')
Tc, yc = ct.step_response(sys, t)
ax.plot ( Tc, yc, linewidth=4, label ='with compensation')
ax.set_xlabel ('Time(s)', fontsize =16)
ax.set_ylabel ('Output', fontsize =16)
plt.show()