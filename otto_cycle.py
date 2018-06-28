#adiabat compression
#isochor 
#adiabat expansion
#isovhor

import matplotlib.pyplot as plt
import numpy as np

v_ratio=15
T_max=1200
T_1=293
gamma=1.4

dv=0.1
#P1V1^(gamma)=PnVn^(gamma)
#T1V1^(gamma-1)=TnVn^(gamma-1)
#T1/P1V1=T2/P2V2
adiabat_constant_1=1
v_1_2=np.linspace(1,1/v_ratio,100)
p_1_2=[adiabat_constant_1/(x**gamma) for x in v_1_2]
plt.plot(v_1_2,p_1_2)
T_2=(T_1*p_1_2[-1]*v_1_2[-1])/(v_1_2[0]*p_1_2[0])


#P2=P1T2/T1
p_3=p_1_2[-1]*T_max/T_2
p_2_3=np.linspace(p_1_2[-1],p_3,100)
v_2_3=[1/v_ratio for x in p_2_3]
plt.plot(v_2_3,p_2_3)


adiabat_constant_2=p_2_3[-1]*v_2_3[-1]**(gamma)
v_3_4=np.linspace(v_2_3[-1],1,100)
p_3_4=[adiabat_constant_2/(x**gamma) for x in v_3_4]

plt.plot(v_3_4,p_3_4)

p_4_1=np.linspace(p_3_4[-1],1,100)
v_4_1 = np.ones(len(p_4_1))
T_4=(T_max*p_3_4[-1]*v_3_4[-1])/(v_3_4[0]*p_3_4[0])
plt.plot(v_4_1,p_4_1)
plt.show()
area_1=np.trapz(p_1_2,v_1_2)+np.trapz(p_2_3,v_2_3)+np.trapz(p_3_4,v_3_4)+np.trapz(p_4_1,v_4_1)


#step 1-2 adiabat comp
s_0=0
t_1_2=np.linspace(T_1,T_2,100)
s_1_2=[s_0 for x in t_1_2]
plt.plot(s_1_2,t_1_2)
print("T1={0:.2f},S1={1:.2f}".format(t_1_2[0],s_1_2[0]))

#step 2-3 isochoric heating
t_2_3=np.linspace(T_2,T_max,len(p_2_3))
s_2=s_1_2[-1]
s_2_3=np.empty(len(p_2_3))
s_2_3[0]=s_2
#Tn=Ti*Pn*Vn/(Pi*Vi)
t_2_3=[(T_2*p_2_3[i]*v_2_3[i])/(v_2_3[0]*p_2_3[0]) for i in range(len(p_2_3))]
for i in range(1,len(p_2_3)):
    s_2_3[i]=s_2_3[i-1]+(5/2)*(p_2_3[i]-p_2_3[i-1])*(v_2_3[i])/(t_2_3[i])

plt.plot(s_2_3,t_2_3)
print("T2={0:.2f},S2={1:.4f}".format(t_2_3[0],s_2_3[0]))

#step 3-4 adiabat exp
s_3=s_2_3[-1]
t_3_4=np.linspace(T_max,T_4,100)
s_3_4=[s_3 for x in t_3_4]
plt.plot(s_3_4,t_3_4)
print("T3={0:.2f},S3={1:.4f}".format(t_2_3[-1],s_2_3[-1]))

#step 4-1 isochoric cooling
t_4_1=np.linspace(T_4,T_1,len(p_4_1))
s_4=s_3_4[-1]
s_4_1=np.empty(len(p_4_1))
s_4_1[0]=s_4
#Tn=Ti*Pn*Vn/(Pi*Vi)
t_4_1=[(T_4*p_4_1[i]*v_4_1[i])/(v_4_1[0]*p_4_1[0]) for i in range(len(p_2_3))]
for i in range(1,len(p_2_3)):
    s_4_1[i]=s_4_1[i-1]+(5/2)*(p_4_1[i]-p_4_1[i-1])*(v_4_1[i])/(t_4_1[i])

plt.plot(s_4_1,t_4_1)
print("T4={0:.2f},S4={1:.4f}".format(t_4_1[0],s_4_1[0]))

#print("{0:.4f} {1:.4f} {2:.4f} {3:.4f}".format(np.trapz(t_1_2,s_1_2),np.trapz(t_2_3,s_2_3),np.trapz(t_3_4,s_3_4),np.trapz(t_4_1,s_4_1)))
area_2=np.trapz(t_2_3,s_2_3)#+np.trapz(t_4_1,s_4_1)
print("area 1={0:.4f}, area 2={1:.4f}".format(area_1,area_2))
print("eff={0:.2f}".format(area_1/area_2))

plt.show()