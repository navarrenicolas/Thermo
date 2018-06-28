import matplotlib.pyplot as plt
import numpy as np
import colorsys
from scipy.integrate import quad

kb=1.38064852e-23
def prob(v,m,t):
    return np.sqrt((m/(2*np.pi*(kb)*t))**3)*4*np.pi*(v**2)*np.exp(-(m*(v**2))/(2*(kb)*t))

velocity = np.arange(0,3500,2)
atomic_mass=1.67377e-27
mass_helium=4.002602*atomic_mass
mass_air=28.97*atomic_mass

temp=293.15
# for temp in range (50,300):
#     HSV = [(temp*-0.280e-2+0.84, 1, 1) for x in velocity]
#     RGB = [colorsys.hsv_to_rgb(x,y,z) for x,y,z in HSV]
#     # probability_helium=[prob(x,mass_helium,temp) for x in velocity]
#     probability_air=[prob(x,mass_air,temp) for x in velocity]
#     #scale=5000
#     # plt.scatter(velocity,probability_helium,c=RGB)
#     plt.scatter(velocity,probability_air,c=RGB,s=0.8)


probability_helium=[prob(x,mass_helium,temp) for x in velocity]
probability_air=[prob(x,mass_air,temp) for x in velocity]
plt.plot(velocity,probability_helium,label='helium')
plt.plot(velocity,probability_air,label='air')
plt.grid(True)
plt.legend()
plt.xlabel("Velocity $m/s$")
plt.ylabel("Probability Density Function")
plt.title("Maxwell-Boltzmann Distribution")
plt.show()

def v_rms(t,m):
    return np.sqrt(3*kb*t/m)
def c(v,gamma):
    return gamma*v/np.sqrt(3)

vrms_air=v_rms(293.15,mass_air)
vrms_helium=v_rms(293.15,mass_helium)
c_air=c(vrms_air,7/5)
c_helium=c(vrms_helium,5/3)

print("v_rms(air)={0:.2f}m/s, vrms(helium)={1:.2f}m/s".format(vrms_air,vrms_helium))
print("c(air)={0:.2f}m/s, c(helium)={1:.2f}m/s".format(c_air,c_helium))

