import numpy as np
import matplotlib.pyplot as plt
from fwmi import Hybrid


switch_lam = 2
lam = 532e-9 * switch_lam
theta_t = 1.5 * np.pi / 180
gamma_m = 1.40e9 / switch_lam  # molecular signal spectral width
gamma_a = 50e6 / switch_lam  # aerosol signal spectral width
fopd = 0.1 * switch_lam
t_ref = 20
t = 20
p = 1
f = 0.1

h1 = Hybrid(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p)
n_nu = 100
theta_d = 0.5 * np.pi / 180
d_nu = np.linspace(0.01e9, 0.1e9, n_nu)

t_m = np.zeros(n_nu)
t_a = np.zeros(n_nu)
sdr = np.zeros(n_nu)
opd = h1.opd_exact(theta_t, h1.n1, h1.d1, h1.n2, h1.d2, h1.n3, h1.d3)
for i in range(n_nu):
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(opd), phase_dev=d_nu[i])
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(opd), phase_dev=d_nu[i])
    sdr[i] = t_m[i] / t_a[i]

fig, ax = plt.subplots()
ax.plot(d_nu / 1e9, sdr, color='black')
ax.grid(True)
ax.set_xlabel(r"Locking Error (GHz)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Locking Error")
plt.show()
