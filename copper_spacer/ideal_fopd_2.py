"""

"""
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from fwmi import Metal
from glass_library import LITHOSIL_Q, N_SF66

mol_fwhm_1064 = 1.4e9
aer_fwhm_1064 = 100e6
switch_lam = 1
lam = 1064e-9 * switch_lam
theta_t = 1.5 * np.pi / 180
gamma_m = mol_fwhm_1064 / 2  # molecular signal spectral half width (1 / e)
gamma_a = aer_fwhm_1064 / 2  # aerosol signal spectral half width (1 / e)
fopd = 0.15 * switch_lam
t_ref = 20
t = 20
p = 1
f = 0.1

# h1 = H1(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p)

n_fopds = 40
theta_d = 0.002
fopds = np.linspace(0.05, 0.3, n_fopds)

t_m = np.zeros(n_fopds)
t_a = np.zeros(n_fopds)
sdr = np.zeros(n_fopds)
fsrs = np.zeros(n_fopds)
theta_t_0 = theta_t
for i in range(n_fopds):
    h1 = Metal(fopds[i], theta_t, gamma_m, gamma_a, lam, t, t_ref, p, d_opd_d_t=lam / 5, glass=LITHOSIL_Q)
    fsrs[i] = h1.fsr(fopds[i])
    t_m[i] = h1.overall_transmittance(theta_d, f, h1.gamma_m, h1.fsr(fopds[i]))
    t_a[i] = h1.overall_transmittance(theta_d, f, h1.gamma_a, h1.fsr(fopds[i]))
    sdr[i] = t_m[i] / t_a[i]


# fig, ax = plt.subplots()
# ax.plot(fsrs / 1e9, sdr, color='black')
# # ax.plot(d_d1 * 1000, fsrs, color='black')
# # ax.set_ylim([0, 400])
# ax.grid(True)
# ax.set_xlabel(r"$\Delta\theta_t$ (degrees)")
# ax.set_ylabel(r"SDR")
# ax.set_title("SDR v. Tilt Angle Variation")
# # ax.text(2.1, 110, r"$\theta_t$($\Delta\theta_t$=0) = {0}$^\circ$".format(round(theta_t_0 * 180 / np.pi, 1)))
# plt.show()
# fig, ax = plt.subplots(3)
# ax[0].plot(fsrs / 1e9, t_m)
# ax[0].set_ylabel(r'T$_m$')
# ax[0].set_title(r'$\gamma_m$ = {0} GHz, $\gamma_a$ = {1} MHz'.format(2 * gamma_m / 1e9, 2 * gamma_a / 1e6))
# ax[1].plot(fsrs / 1e9, t_a)
# ax[1].set_ylabel(r'T$_a$')
# ax[2].plot(fsrs / 1e9, sdr)
# ax[2].set_xlabel('FSR (GHz)')
# ax[2].set_ylabel(r'T$_m$ / T$_a$')
# plt.show()
fig, ax = plt.subplots(2)
ax[0].plot(fsrs / 1e9, t_m)
ax[0].grid(True)
ax[0].set_ylabel(r'T$_m$')
ax[0].set_title(r'$\gamma_m$ = {0} GHz, $\gamma_a$ = {1} MHz'.format(2 * gamma_m / 1e9, 2 * gamma_a / 1e6))
ax[1].plot(fsrs / 1e9, sdr)
ax[1].grid(True)
ax[1].set_xlabel('FSR (GHz)')
ax[1].set_ylabel(r'T$_m$ / T$_a$')
plt.show()

