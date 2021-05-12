import numpy as np
import matplotlib.pyplot as plt
from fwmi import Hybrid
from glass_library import *

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

h1 = Hybrid(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p, d_opd_d_t=lam / 5, pure_arm=LITHOSIL_Q, hybrid_arm=BK7)
n_angles = 30
theta_d = np.linspace(0.1, 3, n_angles) * np.pi / 180

t_m = np.zeros(n_angles)
t_a = np.zeros(n_angles)
sdr = np.zeros(n_angles)
for i in range(n_angles):
    t_m[i] = h1.overall_transmittance(theta_d[i], f, h1.gamma_m, h1.fsr(fopd))
    t_a[i] = h1.overall_transmittance(theta_d[i], f, h1.gamma_a, h1.fsr(fopd))
    sdr[i] = t_m[i] / t_a[i]

fig, ax = plt.subplots()
ax.plot(theta_d * 1000, sdr, color='black')
ax.grid(True)
ax.set_xlabel(r"$\theta_d$ (mrad)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Half Divergent Angle")
plt.show()
