"""
Plot the spectral discrimination ratio versus the aerosol backscatter spectral width.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015).

All units are SI.
"""
import numpy as np
import matplotlib.pyplot as plt
import design_library as dlib

gamma_m = 1.40e9 / 2  # molecular signal spectral width
fopd = 0.2

gamma_a = np.linspace(15, 50, 201) * 1e6
t_m = dlib.transmittance_simple(gamma_m, fopd)
t_a = dlib.transmittance_simple(gamma_a, fopd)
sdr = dlib.sdr_simple(gamma_m, gamma_a, fopd)

fig, ax = plt.subplots()
ax.plot(gamma_a / 1e6, sdr)
ax.set_xlabel("Aerosol Linewidth (MHz)")
ax.set_ylabel("SDR")
ax.grid(True)
plt.text(35, 725, r"$\gamma_m$ = {0} Ghz, FOPD = {1} m".format(round(gamma_m / 1e9, 2), fopd))
# plt.show()
plt.savefig("sdr_v_aerosol_width_1064.pdf")
plt.close()