import numpy as np
import matplotlib.pyplot as plt
from fwmi import Metal
from glass_library import LITHOSIL_Q, N_SF66


def hm(dist):
    """

    :param dist:
    :return:
    """
    max = np.amax(dist)
    half_max = max / 2
    max_ind = np.argmax(dist)
    return_index = None
    return_half_max = None
    for index, value in enumerate(dist):
        if value >= half_max:
            return_index = index
            return_half_max = value
            break
    return return_index, return_half_max

switch_lam = 2
lam = 532e-9 * switch_lam
theta_t = 1.5 * np.pi / 180
gamma_m = 1.40e9 / switch_lam  # molecular signal spectral width
gamma_a = 100e6 / switch_lam  # aerosol signal spectral width
fopd = 3e8 / 2e9 # * switch_lam
t_ref = 20
t = 20
p = 1
f = 0.1

n_nu = 100
n_fopds = 40
theta_d = 0.002
fopds = np.linspace(0.05, 0.3, n_fopds)
d_nu = np.linspace(-0.1e9, 0.1e9, n_nu)
half_maxs = np.zeros(n_fopds)
hwhms = np.zeros(n_fopds)
fsrs = np.zeros(n_fopds)

t_m = np.zeros(n_nu)
t_a = np.zeros(n_nu)
sdr = np.zeros(n_nu)
for j in range(n_fopds):
    h1 = Metal(fopds[j], theta_t, gamma_m, gamma_a, lam, t, t_ref, p, d_opd_d_t=lam / 5, glass=N_SF66, metal_a=None)
    fsrs[j] = h1.fsr(fopds[j])
    for i in range(n_nu):
        t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(fopds[j]), phase_dev=d_nu[i])
        t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(fopds[j]), phase_dev=d_nu[i])
        sdr[i] = t_m[i] / t_a[i]
    info = hm(sdr)
    half_maxs[j] = info[1]
    hwhms[j] = d_nu[info[0]]
    print(info, fsrs[j], d_nu[info[0]] / 1e6)

fig, ax = plt.subplots()
ax.plot(fsrs / 1e9, half_maxs, color='black')
ax.grid(True)
ax.set_xlabel(r"FSR (GHz)")
ax.set_ylabel(r"Half Max")
ax.set_title("Locking Error Half Max v FSR")
ax.text(3, 45, "Half Width = {0} MHz".format(round(abs(hwhms[0]) / 1e6, 1)))
plt.show()
