"""
Determine the ideal fixed optical path difference for a FWMI, with no imperfections or axis tilt.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015) and corresponds to
Fig. 2 in that paper.

All units are SI.
"""
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import design_library as dlib

gamma_m = 1.40e9 / 2  # molecular signal spectral width
gamma_a = 100e6 / 2  # aerosol signal spectral width
gamma_l = 0e6  # laser spectral width

fopd = np.linspace(0.05, 0.3, 201)
t_m = dlib.transmittance_simple(gamma_m, fopd)
t_a = dlib.transmittance_simple(gamma_a + gamma_l, fopd)
sdr = dlib.sdr_simple(gamma_m, gamma_a + gamma_l, fopd)

specific_fopd = 0.1
t_m_100mm = dlib.transmittance_simple(gamma_m, specific_fopd)
sdr_100mm = dlib.sdr_simple(gamma_m, gamma_a + gamma_l, specific_fopd)
print(t_m_100mm, sdr_100mm)

fig, ax = plt.subplots(2)
ax[0].plot(const.c / fopd / 1e9, t_m)
ax[0].set_ylabel(r'T$_m$')
ax[0].set_title(r'$\gamma_m$ = {0} GHz, $\gamma_a$ = {1} MHz'.format(2 * gamma_m / 1e9, 2 * gamma_a / 1e6))
ax[1].plot(const.c / fopd / 1e9, sdr)
ax[1].set_xlabel('FSR (GHz)')
ax[1].set_ylabel(r'T$_m$ / T$_a$')
plt.show()