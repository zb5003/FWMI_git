"""
Determine the thicknesses of the three materials in the arms of the FWMI.
Follows [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015).
The material indices are as follows:
 -1: pure arm glass
 -2: hybrid arm glass
 -3: hybrid arm air.

All units are SI and radians unless otherwise stated.
"""
import numpy as np
import numpy.linalg as la
from numpy import sin, sqrt, asarray
import design_library as dlib


def n_air(lam, t, p):
    """
    Calculate the index of refraction of air.

    :param lam: Float. The wavelength at which to calculate the refractive index, in microns.
    :param t: Float. The temperature, in celsius, at which to calculate the refractive index.
    :param p: Float. The pressure, in atmospheres, at which to calculate the refractive index.
    :return: Float. The refractive index of air.
    """
    return 1 + (6432.8 + 2949810 * lam**2 / (146 * lam**2 - 1) + 25540 * lam**2 / (41 * lam**2 - 1)) / (1 + (t - 15) * 3.4785e-3) * 1e-8 * p


def beta_air(lam, t, p):
    """
    Calculate the temperature derivative of the refractive index of air.

    :param lam: Float. The wavelength at which to calculate the refractive index, in microns.
    :param t: Float. The temperature, in celsius, at which to calculate the refractive index.
    :param p: Float. The pressure, in atmospheres, at which to calculate the refractive index.
    :return: Float. dn_air / dT
    """
    return -3.4785e-3 / (1 + (t - 15) * 3.4785e-3) * (n_air(lam, t, p) - 1)


def delta_n(n_ref, d0, d1, d2, e0, e1, lam_tk, lam, t, t_ref=20):
    """
    Calculate the change in refractive index, with respect to the reference index, as a function of temperature
    deviation from the reference temperature.

    This equation can be found in [1] as well as chapter 22 of the Zemax user manual.
    :param n_ref: Float. The reference refractive index at the reference temperature at the wavelength of interest.
    :param d0: Float.
    :param d1: Float.
    :param d2: Float.
    :param e0: Float.
    :param e1: Float.
    :param lam_tk: Float. Reference wavelength, in microns. Determined at the reference temperature, but is a fit
                   parameter and not an explicit function of index or temperature.
    :param lam: Float. Wavelength of interest, in microns.
    :param t: Float. The temperature at which to calculate the refractive index change.
    :param t_ref. Float. The temperature at which the reference refractive index was measured. Default is 20 Celsius.
    :return: Float. The temperature derivative of the refractive index.
    """
    dt = t - t_ref
    return (n_ref**2 - 1) / (2 * n_ref) * (d0 * dt + d1 * dt**2 + d2 * dt**3 + (e0 * dt + e1 * dt**2) / (lam**2 - lam_tk**2))


def beta(n_ref, d0, d1, d2, e0, e1, lam_tk, lam, t, t_ref=20):
    """
    Calculate the temperature derivative of the refractive index of glass.

    :param n_ref: Float. The reference refractive index at the reference temperature at the wavelength of interest.
    :param d0: Float.
    :param d1: Float.
    :param d2: Float.
    :param e0: Float.
    :param e1: Float.
    :param lam_tk: Float. Reference wavelength, in microns. Determined at the reference temperature, but is a fit
                   parameter and not an explicit function of index or temperature.
    :param lam: Float. Wavelength of interest, in microns.
    :param t: Float. The temperature at which to calculate the refractive index change.
    :param t_ref. Float. The temperature at which the reference refractive index was measured. Default is 20 Celsius.
    :return: Float. The temperature derivative of the refractive index.
    """
    dt = t - t_ref
    return (n_ref**2 - 1) / (2 * n_ref) * (d0 + 2 * d1 * dt + 3 * d2 * dt**2 + (e0 + 2 * e1 * dt) / (lam**2 - lam_tk**2))


def sellmeier_1(lam, k1, l1, k2, l2, k3, l3):
    """
    Calculate the index of refraction using the Sellmeier 1 dispersion formula in chapter 21 of the Zemax manual.

    These coefficients are measured at a particular reference temperature. To find the index at the wavelength of
    interest at a different temperture, use this formula in conjucntion with delta_n().
    :param lam: Float. The wavelength of interest.
    :param k1: Float.
    :param l1: Float.
    :param k2: Float.
    :param l2: Float.
    :param k3: Float.
    :param l3: Float.
    :return: Float. The index of refraction at the wavelength of interest.
    """
    return sqrt(1 + k1 * lam**2 / (lam**2 - l1) + k2 * lam**2 / (lam**2 - l2) + k3 * lam**2 / (lam**2 - l3))

theta_t = 1.5 * np.pi / 180
fopd = 0.1
lam_of_interest = 532e-3
t = 19  # temperature in Celsius
p = 1  # atmospheric pressure, in atm

# N-SF66
t_ref_nsf66 = 20
n_ref_nsf66 = sellmeier_1(lam_of_interest, 2.024597, 0.0147053225, 0.470187196, 0.0692998276, 2.59970433, 161.817601)  # 0.39 - 2.5 microns
n1 = n_ref_nsf66 + delta_n(n_ref_nsf66, -4.3e-6, 1.15e-8, 4.31e-11, 9.62e-7, 1.62e-9, 3.22e-1, lam_of_interest, t, t_ref_nsf66)
a1 = 5.9e-6
b1 = beta(n_ref_nsf66, -4.3e-6, 1.15e-8, 4.31e-11, 9.62e-7, 1.62e-9, 322e-3, lam_of_interest, t)

# P-SF68
t_ref_psf68 = 20
n_ref_psf68 = sellmeier_1(lam_of_interest, 2.3330067, 0.0168838419, 0.452961396, 0.0716086325, 1.25172339, 118.707479)  # 0.42 - 2.5 microns
n2 = n_ref_psf68 + delta_n(n_ref_psf68, 1.55e-5, 2.3e-8, -3.46e-11, 2.76e-6, 2.93e-9, 297e-3, lam_of_interest, t, t_ref_psf68)
a2 = 8.43e-6  # TCE -30 to 70 C
b2 = beta(n_ref_psf68, 1.55e-5, 2.3e-8, -3.46e-11, 2.76e-6, 2.93e-9, 297e-3, lam_of_interest, t)

n3 = n_air(lam_of_interest, t, p)
a3 = 0
b3 = beta_air(lam_of_interest, t, p)

root1 = sqrt(n1**2 - sin(theta_t)**2)
root2 = sqrt(n2**2 - sin(theta_t)**2)
root3 = sqrt(n3**2 - sin(theta_t)**2)

m1 = 2 * asarray([[n1 * sqrt(1 - (sin(theta_t) / n1)**2), -n2 * sqrt(1 - (sin(theta_t) / n2)**2), -n3 * sqrt(1 - (sin(theta_t) / n3)**2)],
                  [-1 / (2 * root1), 1 / (2 * root2), 1 / (2 * root3)],
                  [a1 * root1 + n1 * b1 / root1, -(a2 * root2 + n2 * b2 / root2), -(a3 * root3 + n3 * b3 / root3)]])
v1 = asarray([fopd, 0, 0])

thickness = la.solve(m1, v1)
print(thickness, [n_ref_nsf66, n_ref_psf68, n3], [n1, n2, n3])
