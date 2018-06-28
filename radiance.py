import numpy as np
import matplotlib.pyplot as plt
from sympy import Eq, Symbol, solve, var
'''
Run the code,
all plots will appear
'''

kb=1.38064852e-23
h=6.626e-34
c=3e8
def formula(lam,T):
    return (2*h*(c**2))/((lam**5)*(np.exp((h*c)/(lam*kb*T))-1))

def formula_ev(e,T):
    return ((e**3))/(((4.1357e-15)**2)*(c**2)*(np.exp((e)/((8.62e-5)*T))-1))

k_b = 1.38 * 10 ** -23  # boltzmann constant
Temp = 5800  # K
h = 6.626 * 10 ** -34  # J s # plank constant
h_ev = 4.135667662 * 10 ** -15  # eV s # plank constant
c = 3 * 10 ** 8  # m/s speed of light
r_earth = 6.371 * 10 ** 6  # m radius of Earth
r_sun = 6.955088 * 10 ** 8  # radius of sun
dis = 1.49597870700 * 10 ** 11  # distance from center to center
Area = np.pi * r_sun ** 2  # area of exposed earth surface
omega = Area / (dis**2)  # stra

spectral_radiance_freq = lambda w: omega * (2 * h * w**3) / (c**2 * (np.exp(h * w / (k_b * Temp)) - 1))

r_e=6.371e6
r_es=149.6e9
nm_ev=6.24150913e18
str_ra=np.pi*(r_e**2)/(r_es**2)
print("Solid angle conversion factor")
print(str_ra)
lam_range=np.linspace(100e-9,3e-6,5000)
irradiance=[formula(x,5800) for x in lam_range]

ev_range=np.arange(0.1,5,0.01)
irradiance_ev=[formula_ev(x,5800)*str_ra for x in ev_range]

irradiance_nm_m_3=[x*omega*1e-9 for x in irradiance]
lam_range_mum=lam_range*(10**6)
plt.plot(lam_range_mum,irradiance_nm_m_3)
plt.title("Spectral Radiance at Earth")
plt.ylabel("$W/m^2 nm$")
plt.xlabel("Wavelength $\mu m$")

plt.show()

energy = np.arange(0.1, 10, 0.1)

plt.plot(energy,spectral_radiance_freq(energy/h_ev)/h_ev)
plt.title("Spectral Radiance at Earth")
plt.ylabel("$W/m^2 ev$")
plt.xlabel("Energy $ev$")

plt.show()
int1=np.trapz(irradiance_nm_m_3,lam_range*(10**9))
print("\nSolar intensity approximation")
print(int1)


def q_3():
    global irradiance_ev,ev_range,lam_range_mum,irradiance_nm_m_3
    lam_range_nm=lam_range_mum*(10**3)
    def a():
        bandgap=1.5
        f_e=np.empty(len(ev_range))
        for i in range(len(ev_range)):
            if ev_range[i]<bandgap:
                f_e[i]=0
            else:
                f_e[i]=irradiance_ev[i]*1.5/ev_range[i]

        plt.plot(ev_range,f_e,label="photvoltaic irradiance")
        plt.plot(ev_range,irradiance_ev,label="solar irradiance")
        plt.title("Solar and Photvoltaic Irradiance for a 1.5ev Bandgap Cell")
        plt.xlabel("Energy $ev$")
        plt.ylabel("$W/m^2 ev$")
        plt.grid(True)
        plt.legend()
        plt.show()
    def b():
        bandgap=1.5
        f_e=np.empty(len(lam_range_nm))
        for i in range(len(lam_range_nm)):
            if 1240/lam_range_nm[i]<bandgap:
                f_e[i]=0
            else:
                f_e[i]=irradiance_nm_m_3[i]*1.5/(1240/lam_range_nm[i])

        plt.plot(lam_range_nm,f_e,label="photvoltaic irradiance")
        plt.plot(lam_range_nm,irradiance_nm_m_3,label="solar irradiance")
        plt.title("Solar and Photvoltaic Irradiance for a 1.5ev Bandgap Cell")
        plt.xlabel("Wavelength $nm$")
        plt.ylabel("$W/m^2 nm$")
        plt.grid(True)
        plt.legend()
        plt.show()
        int1=np.trapz(irradiance_nm_m_3,lam_range_nm)
        int2=np.trapz(f_e,lam_range_nm)
        eff=int2/int1
        print("\nBandgap Efficiency")
        print(eff)
    a()
    b()
        

def q_4():
    r=0.05

    kc=5
    ep=1
    solar_intensity=1000 #W/m^2
    t_amb=293.15
    s_a=4*np.pi*(r**2)
    kb=1.38064852e-23
    var('t')
    eqn = Eq(1000*np.pi*(r**2),s_a*5*(t-293.15)+s_a*((t**4)-(t_amb**4))*kb)
    print("\nBall temp solutions, take the reasonable one:")
    print(solve(eqn))

q_3()
q_4()