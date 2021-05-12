"""
The width and exact shape depend on the original theta_t because the glass dimensions are
determined by this theta_t.
"""
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

n_theta = 400
theta_d = 0.002
d_theta_t = np.linspace(-6, 6, n_theta) * np.pi / 180

t_m = np.zeros(n_theta)
t_a = np.zeros(n_theta)
sdr = np.zeros(n_theta)
fsrs = np.zeros(n_theta)
theta_t_0 = theta_t
for i in range(n_theta):
    h1.theta_t = theta_t_0 + d_theta_t[i]
    opd = h1.opd_exact(theta_t,
                   h1.n1,
                   h1.d1,
                   h1.n2,
                   h1.d2,
                   h1.n3,
                   h1.d3)
    fsrs[i] = h1.fsr(opd)
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(opd))
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(opd))
    sdr[i] = t_m[i] / t_a[i]


fig, ax = plt.subplots()
ax.plot(d_theta_t * 1000, sdr, color='black')
# ax.plot(d_d1 * 1000, fsrs, color='black')
# ax.set_ylim([0, 400])
ax.grid(True)
ax.set_xlabel(r"$\Delta\theta_t$ (mrad)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Tilt Angle Variation")
ax.text(30, 110, r"$\theta_t$($\Delta\theta_t$=0) = {0} mrad".format(round(theta_t_0 * 1000, 1)))
plt.show()
