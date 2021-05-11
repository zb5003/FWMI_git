"""
Calculate the optical path difference variation for a given tilt angle versus incident angle, theta.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015), Eq. 2 and corresponds to
Fig. 3a and 3c in that paper.

All units are SI.
"""
import numpy as np
import matplotlib.pyplot as plt
from fwmi import H1


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

h1 = H1(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p)

theta_t = 1.5 * np.pi / 180
theta = np.linspace(-5, 5, 25) * np.pi / 180

opd = h1.opd_exact(theta, h1.n1, h1.d1, h1.n2, h1.d2, h1.n3, h1.d3)
opd_variation = (opd - fopd) / lam

fig, ax = plt.subplots()
ax.plot(theta * 180 / np.pi, opd_variation, color='black')
ax.grid(True)
ax.set_xlabel(r"$\theta$ (Degrees)")
ax.set_ylabel(r"$\Delta$OPD ($\lambda$)")
ax.set_title("OPD Variation v. Incident Angle")
plt.show()