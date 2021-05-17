"""

"""
import numpy as np
import matplotlib.pyplot as plt
from fwmi import Metal
from glass_library import LITHOSIL_Q, N_SF66


switch_lam = 2
lam = 532e-9 * switch_lam
theta_t = 1.5 * np.pi / 180
gamma_m = 1.40e9 / switch_lam  # molecular signal spectral width
gamma_a = 50e6 / switch_lam  # aerosol signal spectral width
fopd = 0.15 * switch_lam
t_ref = 20
t = 20
p = 1
f = 0.1

h1 = Metal(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p, d_opd_d_t=lam / 5, glass=N_SF66, metal_a=17e-6)
print(h1.metal_a, h1.n_glass)

n_T = 1000
theta_d = 0.02
d_T = np.linspace(-2, 2, n_T)
print(h1.d_air, h1.d_glass)

t_m = np.zeros(n_T)
t_a = np.zeros(n_T)
sdr = np.zeros(n_T)
fsrs = np.zeros(n_T)
opds = np.zeros(n_T)
theta_t_0 = theta_t
for i in range(n_T):
    n1 = h1.generate_n_glass(h1.t + d_T[i])
    d1 = h1.d_glass_thermal_expansion(h1.t + d_T[i])
    n2 = h1.generate_n_air(h1.t + d_T[i])
    d2 = h1.d_air_thermal_expansion(h1.t + d_T[i])
    opd = h1.opd_exact_pure(theta_t, n1, d1, n2, d2)
    opds[i] = opd
    fsrs[i] = h1.fsr(opd)
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(opd))
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(opd))
    sdr[i] = t_m[i] / t_a[i]


fig, ax = plt.subplots()
# ax.plot(d_T, t_a, color='black')
ax.plot(d_T, sdr, color='black')
# ax.plot(d_T, (opds - fopd) / lam, color='black')
# ax.plot(d_d1 * 1000, fsrs, color='black')
# ax.set_ylim([0, 400])
ax.grid(True)
ax.set_xlabel(r"Temperature Variation ($^\circ$C)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Temperature Variation")
ax.text(2.1, 110, r"$\theta_t$($\Delta\theta_t$=0) = {0}$^\circ$".format(round(theta_t_0 * 180 / np.pi, 1)))
plt.show()
