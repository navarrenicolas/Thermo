import numpy as np
import matplotlib.pyplot as plt
import colorsys
'''
1. Using the barymetric equation and the dry adiabatic lapse rate, 𝛤𝑑, make a simple finite difference
calculation to find the rate of the change of pressure and density with altitude. Start at 15 C, 101 kPa.
Plot a graph of pressure vs. altitude, and read off the scale height (where the pressure is 1/𝑒 of the sea
level value)?

p**(1/gamma-1)*T=c
'''

#Run this code, the first plot is pressure vs altitude, once you close it  second one will apppear, that one is density vs altitude. 
P0=101
PMin=P0*(1/np.e)
Gamma=1.4
R=8.314
m=28.9
g=9.81
T0=15+273.15
alt_max=40
dx=0.001
steps=int(alt_max/dx)

alt=np.arange(0,alt_max,dx)

dz=alt[1]-alt[0]
temps=np.zeros(steps)
temps[0]=T0
pt_const=(P0**(1/Gamma-1))*T0
dT_dZ=(-m*g/R)*(1-1/Gamma)


for i in range(1,len(temps)):
    temps[i]=(temps[i-1])+(dT_dZ)*dz

pressure=[(pt_const/x)**(1/(1/Gamma-1)) for x in temps if (pt_const/x)**(1/(1/Gamma-1))>=PMin]
print(pt_const)
alt=alt[0:len(pressure)]
temps=temps[0:len(alt)]
mid=len(pressure)//2
##0=m(T[0])+b
##pi=m(T[-1])+b
HSV = [((-0.00806)*x+2.3228, 1, 1) for x in temps]
RGB = [colorsys.hsv_to_rgb(x,y,z) for x,y,z in HSV]
plt.scatter(alt,pressure,s=1,c=RGB)
plt.plot(alt[0],pressure[0],'-ro',c=RGB[0])
plt.plot(alt[-1],pressure[-1],'-ro',c=RGB[-1])
plt.plot(alt[mid],pressure[mid],'-ro',c=RGB[mid])
plt.text(alt[0]+0.3,pressure[0],'T={0:.2f}K\nz={1:.2f}km'.format(temps[0],alt[0]))
plt.text(alt[mid]+0.3,pressure[mid],'T={0:.2f}K\nz={1:.2f}km'.format(temps[mid],temps[mid]))
plt.text(alt[-1]+0.3,pressure[-1],'T={0:.2f}K\nz={1:.2f}km'.format(temps[-1],alt[-1]))
plt.axis([-0.5,10,30,110])
plt.legend()
plt.title('Pressure vs Altitude and Tempurature in Earth Atmosphere')
plt.xlabel('Altitude $(km)$')
plt.ylabel('Pressure $(kPa)$')
plt.grid(True)
plt.show()

cg10 = lambda z: -m ** 2 * P0 * T0 ** ((-1) * 0.14e1 / (0.14e1 - 1)) * (-(0.14e1 - 1) * z * m * g + R * 0.14e1 * T0) ** ((2 - 0.14e1) / (0.14e1 - 1)) * g * R ** ((-1) * 0.14e1 / (0.14e1 - 1)) * 0.14e1 ** (-0.1e1 / (0.14e1 - 1))
cg11 = lambda z: m * P0 * T0 ** ((-1) * 0.14e1 / (0.14e1 - 1)) * (-(0.14e1 - 1) * z * m * g + R * 0.14e1 * T0) ** (0.1e1 / (0.14e1 - 1)) * R ** ((-1) * 0.14e1 / (0.14e1 - 1)) * 0.14e1 ** (-0.1e1 / (0.14e1 - 1))


#plt.scatter(alt,[cg10(x) for x in alt],label='$p\'(z)$',s=pressure,c=RGB)
density=[cg11(x) for x in alt]
plt.scatter(alt,density,s=[x*10 for x in pressure],c=RGB)
plt.plot(alt,density,label='$p(z)$',color='black')
plt.plot(alt[0],density[0],'-ro')
plt.plot(alt[-1],density[-1],'-ro')
plt.plot(alt[mid],density[mid],'-ro')
plt.axis([-0.5,9.75,0.5,1.4])
plt.text(alt[0]+0.3,density[0],'T={0:.2f}K\nz={1:.2f}km\nP={2:.3f}kPa\np={3:.2f}$kg/m^3$'.format(temps[0],alt[0],pressure[0],density[0]),bbox=dict(facecolor='white', alpha=0.5))
plt.text(alt[mid]+0.3,density[mid],'T={0:.2f}K\nz={1:.2f}km\nP={2:.3f}kPa\np={3:.2f}$kg/m^3$'.format(temps[mid],temps[mid],pressure[mid],density[0]),bbox=dict(facecolor='white', alpha=0.5))
plt.text(alt[-1]+0.3,density[-1],'T={0:.2f}K\nz={1:.2f}km\nP={2:.3f}kPa\np={3:.2f}$kg/m^3$'.format(temps[-1],alt[-1],pressure[-1],density[-1]),bbox=dict(facecolor='white', alpha=0.5))
plt.xlabel('Altitude $(km)$')
plt.ylabel('Density $(kg/m^3)$')
plt.title("Air Density with Respect to Altitude, Tempurature, and Pressure")
plt.grid(True)
plt.legend()
plt.show()