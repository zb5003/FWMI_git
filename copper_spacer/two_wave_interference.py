import numpy as np
import matplotlib.pyplot as plt
from fwmi import Metal
from glass_library import LITHOSIL_Q, N_SF66


def snell_sin_to_cos(n, theta):
    return np.sqrt(1 - np.sin(theta) ** 2 / n ** 2)


def tilt_opd(n1, n2, d1, d2, fopd, dt, theta_0=0, lm=0.02, alpha=15e-6):
    """

    :param n1:
    :param n2:
    :param d1:
    :param d2:
    :param fopd:
    :param dt:
    :param theta_0:
    :param lm:
    :param alpha:
    :return:
    """
    path1 = 2 * n1 * d1 * snell_sin_to_cos(n1, theta_0)
    d_tilt = d2**2 * alpha * dt / lm
    ds = fopd / 2
    path2_coef = np.sqrt(ds**2 + d_tilt**2)
    path2 = 2 * n2 * path2_coef * np.cos(theta_0 + np.arctan(d_tilt / ds))
    return path2  #path1 - path2


def tilt_opd2(n1, n2, d1, d2, theta_0, theta_tilt):
    """

    :param n1:
    :param n2:
    :param d1:
    :param d2:
    :param theta_0:
    :param theta_tilt:
    :return:
    """
    term1 = 2 * d1 * n1 * snell_sin_to_cos(n1, theta_0)
    term2 = d2 * n2 * snell_sin_to_cos(n2, theta_0)
    theta2 = np.arcsin(np.sin(theta_0) / n2)
    term3 = d2 * n2 * snell_sin_to_cos(n2, theta2 + theta_tilt)
    return term1 - term2 - term3


lam = 1064e-9
theta_t = 1.5 * np.pi / 180
gamma_m = 1.40e9 / 2  # molecular signal spectral width
gamma_a = 100e6 / 2  # aerosol signal spectral width
fopd = 3e8 / 2e9
t_ref = 20
t = 20
p = 1
f = 0.1

h1 = Metal(fopd, theta_t, gamma_m, gamma_a, lam, t, t_ref, p, d_opd_d_t=lam / 5, glass=N_SF66, metal_a=None)
thetas = np.linspace(-0.2, 0.2, 1001) + 0.025
d_temp = np.linspace(0, 1, 101)
opds = h1.opd_exact_pure(thetas, h1.n_glass, h1.d_glass, h1.n_air, h1.d_air)
# opds2 = tilt_opd(h1.n_glass, h1.n_air, h1.d_glass, h1.d_air, fopd, dt=d_temp)
opds3 = tilt_opd(h1.n_air, h1.n_air, h1.n_glass * h1.d_glass, h1.d_air, fopd, dt=0.1, theta_0=thetas)
opds4 = tilt_opd2(h1.n_glass, h1.n_air, h1.d_glass, h1.d_air, thetas, np.arctan(h1.d_air * 15e-6 * 1 / 0.02))
transfer = 1 + np.cos(2 * np.pi / lam * opds)
print('asdfasd', h1.d_air * 15e-6 * 0.1 / 0.02)
# transfer2 = 1 + np.cos(2 * np.pi / lam * opds2)
transfer3 = 1 + np.cos(2 * np.pi / lam * opds3)
transfer4 = 1 + np.cos(2 * np.pi / lam * opds4)

# fig, ax = plt.subplots(2)
# ax[0].plot(thetas * 1000, transfer4)
# ax[1].plot(thetas * 1000, transfer)
# ax[1].set_xlabel('Incident Angle (mrad)')
# ax[0].set_title('No Tilt')
# plt.show()

fig, ax = plt.subplots()
ax.plot(thetas * 1000, transfer, color='red', linestyle='dashed')
ax.plot(thetas * 1000, transfer4, color='black')
ax.set_xlabel('Incident Angle (mrad)')
ax.set_title('No Tilt')
plt.show()
