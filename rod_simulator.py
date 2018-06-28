import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


#Constants
L=0.50 #length of rod
T=3000 #time of simulation
Nx=50 #spacial divisions
Nt=2000 #time divisons
k=16.2
c=502
rho=8000
radius=0.005
t_ambient=20
avgerage=False #Typo on average but set this to true to see the average temp of rod vs time (only if animate is false) 
boltz=5.67e-8
kel=273.15

#heat loss constants
convection_const=30.3
emissivity=0.4

#initialize arrays, time array is actually useless
x = np.linspace(0, L, Nx+1)    
dx = x[1] - x[0]

t = np.linspace(0, T, Nt+1)    
dt = t[1]-t[0]

avg_t=np.zeros(Nt+1)

#more constants/initial conditions
alpha = (k*dt)/((dx**2)*c*rho)
power=2
power_const =(power*dt)/(c*rho*np.pi*(radius**2)*dx)
u_1 = 20*np.ones(Nx+1)

convection_loss_constant=(2*convection_const*dt)/(c*rho*radius)
radiation_loss_constant=(2*dt*emissivity*boltz)/(c*rho*radius)

convection_end_loss_constant=(dt*convection_const)/(c*rho*dx) 
radiation_end_loss_constant=(dt*boltz*emissivity)/(c*rho*dx) 

#latches, time keeping, and indexing constants 
done=False
et=0

def apply_alogorithm(temp_profile):
    temp_profile_copy=np.copy(temp_profile)
    #Finite difference algorithm for loss

    #inside conduction losses
    temp_profile[1:-1] +=alpha*(temp_profile_copy[0:-2]-2*temp_profile_copy[1:-1]+temp_profile_copy[2:])
    
    #conduction losses ends
    temp_profile[0]-=alpha*(temp_profile_copy[0]-temp_profile_copy[1])
    temp_profile[-1]+=alpha*(temp_profile_copy[-2]-temp_profile_copy[-1])

    #power gain
    temp_profile[0]+=power_const
    
    #convection loss
    temp_profile[0]-=convection_end_loss_constant*(temp_profile_copy[0]-t_ambient)
    temp_profile[:]-=convection_loss_constant*(temp_profile_copy[:]-t_ambient)
    temp_profile[-1]-=convection_end_loss_constant*(temp_profile_copy[-1]-t_ambient)

    #radiation_loss
    temp_profile[0]-=radiation_end_loss_constant*(((temp_profile_copy[0]+kel)**4)-((t_ambient+kel)**4))
    temp_profile[:]-=radiation_loss_constant*(((temp_profile_copy[:]+kel)**4)-((t_ambient+kel)**4))
    temp_profile[-1]-=radiation_end_loss_constant*(((temp_profile_copy[-1]+kel)**4)-((t_ambient+kel)**4))
    
    return temp_profile

avg_t[0]=sum(u_1)/len(u_1)
for i in range(1,Nt+1):
    #Finite difference algorithm
    u_1=apply_alogorithm(u_1)
    et+=dt
    avg_t[i]=sum(u_1)/len(u_1)
    if((avg_t[i]-avg_t[i-1])/avg_t[i]<0.0000001):
        break
    
        
if not(avgerage):
    plt.title('Steady State Temperature in an Aluminum Rod\nheated by {0:.0f}$W$ of power, in a {1:.0f}$^o C$ Environment'.format(power,t_ambient))
    plt.plot(x,u_1,label='Temperature profile')
    plt.xlabel('Position m')
    plt.ylabel('Temperature $^oC$')
    plt.text(0.01,u_1[0]-1, 'Average Temperature = {0:.2f}$^o C$'.format(sum(u_1)/len(u_1)))
    plt.text(0.01,u_1[0]-3, 'Time elapsed = {0:.2f}s'.format(et))
    
    #convection power loss
    #2 pi r  kc dx (T-Tamb)
    conv_loss=sum([2*np.pi*radius*dx*convection_const*(temp-t_ambient) for temp in u_1])
    conv_loss_ends=np.pi*convection_const*(radius**2)*((u_1[0]-t_ambient)+(u_1[-1]-t_ambient))

    radiation_loss=sum([2*np.pi*radius*emissivity*boltz*dx*(((temp+kel)**4)-((t_ambient+kel)**4)) for temp in u_1])
    radiation_loss_ends=np.pi*(radius**2)*emissivity*boltz*\
        ((((u_1[0]+kel)**4)-(t_ambient+kel)**4)\
        +(((u_1[-1]+kel)**4)-((t_ambient+kel)**4)))

    print("{0:.2f}".format(conv_loss+radiation_loss+radiation_loss_ends+conv_loss_ends))

else:
    plt.plot(t[0:int(et/dt)],avg_t[0:int(et/dt)], label='Average Temperature')
    plt.title('Average Temperature vs Time\nof a Heated Aluminum Rod')
    plt.xlabel('Time s')
    plt.ylabel('Temperature $^oC$')
    plt.grid(True)
plt.legend()
plt.grid(True)

plt.show()
