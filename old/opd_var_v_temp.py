"""
Calculate the optical path difference variation for a given tilt angle versus temperature.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015), Eqs. 2 and 18 and corresponds to
Fig. 3b and 3d in that paper.

All units are SI.
"""
import numpy as np
import matplotlib.pyplot as plt
import design_library as dlib




# 532
lam = 532e-9
# from material_thickness.py, t=20
# d1 = 0.05322096
# d2 = 0.01679775
# d3 = 0.01915868
n1_ref = 1.9374328041896338
n2_ref = 2.0209909911203465
n3_ref = 1.0002734414349423
# from material_thickness.py, t=19
d1 = 0.05320249
d2 = 0.01678081
d3 = 0.01915746
# n1_ref = 1.9374320625464374
# n2_ref = 2.0209683813019397
# n3_ref = 1.0002743795480677
# Values from paper
# d1 = 0.053220
# d2 = 0.0167970
# d3 = 0.0191608
# n1_ref = 1.9374328041896338
# n2_ref = 2.0209909911203465
# n3_ref = 1.0002734414349423

# 1064
# lam = 1064e-9
# d1 = 0.11111812
# d2 = 0.03433597
# d3 = 0.04158935
# n1_ref = 1.879832806221337
# n2_ref = 1.959546445859997
# n3_ref = 1.0002692825303892

a1 = 5.9e-6
a2 = 8.43e-6
t = 19
t_ref = 20

theta_t = 1.5 * np.pi / 180
dt = np.linspace(-0.5, 0.5, 25)
n1_0 = n1_ref + dlib.delta_n(n1_ref, -4.3e-6, 1.15e-8, 4.31e-11, 9.62e-7, 1.62e-9, 3.22e-1, lam * 1e6, t)
n1 = n1_ref + dlib.delta_n(n1_ref, -4.3e-6, 1.15e-8, 4.31e-11, 9.62e-7, 1.62e-9, 3.22e-1, lam * 1e6, dt + t)
d1_var = d1 * dlib.thermal_expansion(a1, t + dt, t_ref)
n2_0 = n2_ref + dlib.delta_n(n2_ref, 1.55e-5, 2.3e-8, -3.46e-11, 2.76e-6, 2.93e-9, 297e-3, lam * 1e6, t)
n2 = n2_ref + dlib.delta_n(n2_ref, 1.55e-5, 2.3e-8, -3.46e-11, 2.76e-6, 2.93e-9, 297e-3, lam * 1e6, dt + t)
d2_var = d2 * dlib.thermal_expansion(a2, t + dt, t_ref)
n3_0 = dlib.n_air(lam * 1e6, 20, 1)
n3 = dlib.n_air(lam * 1e6, dt + 20, 1)

fopd = dlib.opd_exact(n1_0, d1, n2_0, d2, n3_0, d3, theta_t)
opd = dlib.opd_exact(n1, d1_var, n2, d2_var, n3, d3, theta_t)
opd_variation = (opd - fopd) / lam
print(fopd)
plt.plot(dt, opd_variation)
plt.show()
