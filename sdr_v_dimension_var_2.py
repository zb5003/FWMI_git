import numpy as np
import matplotlib.pyplot as plt
from fwmi import H1


switch_lam = 1
lam = 532e-9 * switch_lam
theta_t = 1.5 * np.pi / 180
gamma_m = 1.40e9 / switch_lam  # molecular signal spectral width
gamma_a = 50e6 / switch_lam  # aerosol signal spectral width
fopd = 0.1 * switch_lam
t_ref = 20
t = 20
p = 1
f = 0.1

h1 = H1(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p)

n_d1 = 100
theta_d = 0.5 * np.pi / 180
d_d1 = np.linspace(-0.03, 0.03, n_d1) / 1000

t_m = np.zeros(n_d1)
t_a = np.zeros(n_d1)
sdr = np.zeros(n_d1)
fsrs = np.zeros(n_d1)
d1 = h1.d1
d3 = h1.d3
for i in range(n_d1):
    h1.d1 = d1 + d_d1[i]
    h1.d3 = d3 - d_d1[i] * h1.n1
    opd = h1.opd_exact(theta_t, h1.n1, h1.d1, h1.n2, h1.d2, h1.n3, h1.d3)
    fsrs[i] = h1.fsr(opd)
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(opd))
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(opd))
    sdr[i] = t_m[i] / t_a[i]


fig, ax = plt.subplots()
ax.plot(d_d1 * 1000, sdr, color='black')
# ax.plot(d_d1 * 1000, fsrs, color='black')
# ax.set_ylim([0, 400])
ax.grid(True)
ax.set_xlabel(r"$\Delta$d$_1$ (mm)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Pure Glass Length Variation")
# ax.text(0.005, 301, r"d$_1$($\Delta$d$_1$=0) = {0} mm".format(round(100 * d1, 2)))
plt.show()
