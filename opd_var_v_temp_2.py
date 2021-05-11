"""
Calculate the optical path difference variation for a given tilt angle versus temperature.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015), Eqs. 2 and 18 and corresponds to
Fig. 3b and 3d in that paper.

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
print(h1.d1, h1.d2, h1.d3)
# dt = np.linspace(-0.5, 0.5, 25)
dt = np.linspace(-10, 10, 25)
opd = h1.opd_exact(theta_t,
                   h1.generate_n1(t + dt),
                   h1.d1_thermal_expansion(t + dt),
                   h1.generate_n2(t + dt),
                   h1.d2_thermal_expansion(t + dt),
                   h1.generate_n3(t + dt),
                   h1.d3)
opd_variation = (opd - fopd) / lam

fig, ax = plt.subplots()
ax.plot(dt, opd_variation, color='black')
ax.grid(True)
ax.set_xlabel(r"$\Delta$T ($^\circ$C)")
ax.set_ylabel(r"$\Delta$OPD ($\lambda$)")
ax.set_title("OPD Variation v. Temperature")
ax.text(0.21, 0.00001, r"T($\Delta$T=0) = 20 $^\circ$C")
plt.show()
