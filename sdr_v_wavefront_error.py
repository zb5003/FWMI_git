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
n_rms = 30
theta_d = 0.5 * np.pi / 180
rms = np.linspace(0, 0.2, n_rms) * lam

t_m = np.zeros(n_rms)
t_a = np.zeros(n_rms)
sdr = np.zeros(n_rms)
for i in range(n_rms):
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(fopd), rms[i])
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(fopd), rms[i])
    sdr[i] = t_m[i] / t_a[i]

fig, ax = plt.subplots()
ax.plot(rms / lam, sdr, color='black')
ax.grid(True)
ax.set_xlabel(r"Wavefront Error ($\lambda$)")
ax.set_ylabel(r"SDR")
ax.set_title("SDR v. Pure Glass Length Variation")
plt.show()
