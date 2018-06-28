#Brayton Cycle
#1-2 Adibiatic compression
#2-3 Isobaric cooling
#3-4 Adibiatic expansion
#4-1 Isobaric heating
import numpy as np
import matplotlib.pyplot as plt

#efficicency theoretical = 1-T1/T2=1-(P1/P2)^(1-1/gamma)

#input constants
P_ratio=20
T_1=298
P_low=100000
V_1=0.8354
T_3=1273
Gamma=1.4
steps=1000
c_p=5/2

def do_computuation(P_ratio,T_1,P_low,V_1,T_3,Gamma,steps,c_p):
    #computed constants
    p_high=P_low*P_ratio
    p_v_ratio_1=(P_low)*(V_1)**Gamma #pv^(gamma)=c
    v_1_2=np.zeros(steps)
    v_1_2[0]=V_1
    p_1_2=np.linspace(P_low,p_high,steps)

    #v(p)=(c/p)**(1/Gamma)
    v_1_2[1:]=(p_v_ratio_1/p_1_2[1:])**(1/Gamma)
    v_2=v_1_2[-1]
    t_2=T_1*(V_1/v_2)**(Gamma-1)
    v_3=T_3*(v_2/t_2)
    p_2_3=np.ones(steps)*p_high
    v_2_3=np.linspace(v_2,v_3,steps)
    p_3_4=np.linspace(p_high,P_low,steps)
    v_3_4=np.zeros(steps)
    v_3_4[0]=v_3
    p_v_ratio_2=(p_high)*(v_3)**(Gamma)

    #v_next=(c/p_next)**(1/Gamma)
    v_3_4[1:]=(p_v_ratio_2/p_3_4[1:])**(1/Gamma)
    v_4=v_3_4[-1]
    t_4=T_3*(v_3/v_4)**(Gamma-1)
    p_4_1=np.ones(steps)*P_low
    v_4_1=np.linspace(v_4,V_1,steps)
    
    work=np.trapz(y=p_2_3, x=v_2_3)+np.trapz(y=p_3_4, x=v_3_4)+np.trapz(y=p_1_2, x=v_1_2)+np.trapz(y=p_4_1, x=v_4_1)
    eff_optimal=(1-(P_low/p_high)**(1-1/Gamma))*100
    delta_u=c_p*p_high*(v_3-v_2)
    work_2_3=np.trapz(y=p_2_3,x=v_2_3)
    heat_in=delta_u+work_2_3
    eff_num=100*work/(heat_in)

    

    return {'v12':v_1_2,'v23':v_2_3,'v34':v_3_4,
        'v41':v_4_1,'p12':p_1_2,'p23':p_2_3,'p34':p_3_4,'p41':p_4_1,
        't1':T_1,'t2':t_2,'t3':T_3,'t4':t_4,'work':work,'efft':eff_optimal,
        'effnum':eff_num,'du':delta_u,'hin':heat_in}
    

def plot_brayton(data):
    # plt.axes(xLim=(-0.3,2.5),yLim=(0,22.5))
    plt.plot(data['v12'],data['p12'],label='1-2 Adibiatic Compression')
    plt.plot(data['v23'],data['p23'],label='2-3 Isobaric Expansion')
    plt.plot(data['v34'],data['p34'],label='3-4 Adibiatic Expansion')
    plt.plot(data['v41'],data['p41'],label='4-1 Isobaric Compression')
    v_1=data['v12'][0];v_2=data['v23'][0];v_3=data['v34'][0];v_4=data['v41'][0]
    p_low=data['p12'][0];p_high=data['p23'][0]
    plt.xlabel("Volume Ratio")
    plt.ylabel("Pressure Ratio")
    plt.title("Brayton Cycle PV Diagram for a Diatomic Gas")
    plt.text(v_1-0.5,p_low-0.6,"T1={0:d}$^oC$".format(data['t1']),backgroundcolor='white',fontsize=8,style='italic')
    plt.text(v_2-0.4,p_high+0.7,"T2={0:.2f}$^oC$".format(data['t2']),backgroundcolor='white',fontsize=8,style='italic')
    plt.text(v_3,p_high+0.7,"T3={0:d}$^oC$".format(data['t3']),backgroundcolor='white',fontsize=8,style='italic')
    plt.text(v_4,p_low+0.6,"T4={0:.2f}$^oC$".format(data['t4']),backgroundcolor='white',fontsize=8,style='italic')
    plt.grid(True)
    plt.plot(v_1,p_low,'ro')
    plt.plot(v_2,p_high,'ro')
    plt.plot(v_3,p_high,'ro')
    plt.plot(v_4,p_low,'ro')
    plt.text(1.2,12.3,
        "Numerical\n -Net Work Done= {0:.2f}\n -Change in Energy (2-3) ={4:.2f}\n -Heat in = {2:.2f}\n -Efficiency = {3:.2f}%\n\nTheoretical\n -Efficiency = {1:.2f}%".format(data['work'],data['efft'],data['hin'],data['effnum'],data['du'])
        ,fontsize=9,backgroundcolor='white',verticalalignment='top', horizontalalignment='left')
    plt.legend()
    plt.show()

x_range=np.linspace(2,400,70)
data_range=[do_computuation(x,T_1,P_low,V_1,T_3,Gamma,steps,c_p)['effnum'] for x in x_range]
plot_brayton(do_computuation(P_ratio,T_1,P_low,V_1,T_3,Gamma,steps,c_p))
plt.scatter(x_range,data_range,label='Numericaly Simulated Data')
plt.title('Efficiency vs Pressure Ration\nof the Brayton Cycle for a Diatomic Gas')

plt.grid(True)
actual=[100*(1-(1/x)**(1-1/Gamma)) for x in x_range]
plt.plot(x_range,actual,label='Analyticaly Calculated Curve',color='salmon')
plt.xlabel('Pressure Ratio')
plt.ylabel('Efficiency %')
plt.legend()
plt.show()





