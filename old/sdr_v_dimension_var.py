"""
Calculate the change in spectral discrimination ratio versus variation in the glass dimensions.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015), Eqs. 2 and 24.

All units are SI.
"""
import numpy as np
import matplotlib.pyplot as plt
import design_library as dlib


# 532
lam = 532e-9
d1 = 0.05322096
d2 = 0.01679775
d3 = 0.01915868
n1 = 1.9374328041896338
n2 = 2.0209909911203465
n3 = 1.0002734414349423
gamma_m = 1.40e9  # molecular signal spectral width
gamma_a = 50e6  # aerosol signal spectral width

# 1064
# lam = 1064e-9
# d1 = 0.11111812
# d2 = 0.03433597
# d3 = 0.04158935
# n1 = 1.879832806221337
# n2 = 1.959546445859997
# n3 = 1.0002692825303892
# gamma_m = 1.40e9 / 2  # molecular signal spectral width
# gamma_a = 50e6 / 2  # aerosol signal spectral width

theta_t = 1.5 * np.pi / 180 * 0
d_d1 = np.linspace(-0.03, 0.03, 100)*2 #/ 1000

fopd = dlib.opd_exact(n1, d1 + d_d1, n2, d2, n3, d3, theta_t)
sdr = dlib.sdr_simple(gamma_m, gamma_a, fopd)

plt.plot(d_d1 + 0.02548, sdr * 322 / 787)
plt.xlim([-0.03, 0.03])
plt.show()
