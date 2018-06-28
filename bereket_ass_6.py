import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integr
import math

integral = integr.quad
x = np.arange(1, 10000, 1)

k_b = 1.38 * 10 ** -23  # boltzmann constant
Temp = 5800  # K
h = 6.626 * 10 ** -34  # J s # plank constant
h_ev = 4.135667662 * 10 ** -15  # eV s # plank constant
c = 3 * 10 ** 8  # m/s speed of light
r_earth = 6.371 * 10 ** 6  # m radius of Earth
r_sun = 6.955088 * 10 ** 8  # radius of sun
dis = 1.49597870700 * 10 ** 11  # distance from center to center
Area = math.pi * r_sun ** 2  # area of exposed earth surface

omega = Area / (dis**2)  # stra

particles_E = lambda E: np.exp(-E / (k_b * Temp))  # most probable number of particles with energy E

photons_E = lambda E: 1 / (particles_E(E) - 1)  # photons

spectral_radiance = lambda l: omega * (2 * h * c**2) / (((l * 10 ** -9)**5) * (np.exp(h * c / (l * 10**-9 * k_b * Temp)) - 1))
spectral_radiance_freq = lambda w: omega * (2 * h * w**3) / (c**2 * (np.exp(h * w / (k_b * Temp)) - 1))

w = 4

radiance, error = integral(spectral_radiance, 0, 8000)
radiance = radiance / 10 ** 9
power, error = integral(lambda x: radiance, 0, 1000)
print("Radiance: {}".format(radiance))

radiance2, error = integral(spectral_radiance_freq, 0, 8000)
radiance2 = radiance
print("Radiance2: {}".format(radiance2))

# q1 A
plt.figure(figsize=(10, 10))
plt.grid('on')
plt.xlabel("lambda [$nm$]", fontsize=20)
plt.ylabel("radiance [$W/(m^2 nm)$]", fontsize=20)
plt.title("radiance as a function of wavelength", fontsize=30)
plt.plot(x, spectral_radiance(x) / 10 ** 9, linewidth=w)
plt.xlim([0, 4000])
plt.ylim(ymin=0)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


def ratio_wavelength(l, gap):
    gap_l = 1240 / gap  # nm
    vals = []

    for x in l:
        if(x <= gap_l):
            vals.append(x / gap_l)
        else:
            vals.append(0)

    return np.array(vals)


def r_w(l, gap):
    gap_l = 1240 / gap  # nm
    if(l <= gap_l):
        return l / gap_l
    else:
        return 0.0


vals = ratio_wavelength(x, 1.1)
plt.plot(x, vals * spectral_radiance(x) / 10 ** 9, linewidth=w)


# q1 B
plt.figure(figsize=(10, 10))
energy = np.arange(0.1, 10, 0.1)

plt.xlabel("energy [$ev$]", fontsize=20)
plt.ylabel("irradiance [$W/(m^2 ev)$]", fontsize=20)
plt.title("irradiance as a function of energy", fontsize=30)
rp = (energy / h_ev)
plt.plot(energy, spectral_radiance_freq(rp) / h_ev, linewidth=w)
plt.grid('on')
plt.ylim(ymin=0)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim([0, 10])

# part 2
# E_photon = lambda l: h_ev * c / (l * 10 ** 9)   # l is in m returns eV


def ratio_energy(E, gap):
    vals = []
    for x in E:
        if(x >= gap):
            vals.append(gap / x)
        else:
            vals.append(0.0)

    return np.array(vals)


def r_e(e, gap):
    if(e >= gap):
        return gap / e
    else:
        return 0.0


#ratio = lambda x, y: x / y if x <= y else 0
vals = ratio_energy(rp * h_ev, 1.1)

plt.plot(energy, vals * spectral_radiance_freq(rp) / h_ev, linewidth=w)


def efficency(gap):
    # energy_ratio = ratio_energy(rp * h_ev, gap)
    # lambda_ratio = ratio_wavelength(x, gap)
    #energy_1 = [1.0 * spectral_radiance_freq(y / h_ev) / h_ev for y in energy]
    #energy_2 = [r_w(v, gap) * spectral_radiance(v) / 10 ** 9 for v in x]

    en_1 = lambda y: r_e(y, gap) * spectral_radiance_freq(y / h_ev) / h_ev
    en_2 = lambda y: r_w(y, gap) * spectral_radiance(y) / 10 ** 9

    #radiance_1 = np.sum(energy_1)
    #radiance_2 = np.sum(energy_2)
    radiance_1, error = integral(en_1, 0.1, 40)
    radiance_2, error = integral(en_2, 0, 10000)
    # integral

    print(radiance_2)
    print(radiance_1)

    return radiance_2 / radiance, radiance_1 / radiance2

e1, e2 = efficency(1.5)

print("Eff 1: {}, Eff 2: {}".format(e1, e2))


plt.show()
