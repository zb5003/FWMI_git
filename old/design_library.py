"""
Library containing the equations necessary to design a FWMI.
Equations come from [1] Z. Cheng, et al. Opt. Exp. 23, 12117 (2015).

All angles are in radians.
"""
import numpy as np
from scipy.integrate import nquad
import scipy.constants as const


def fsr(opd):
    """
    Calculate the free spectral range of the MI.

    :param opd: Float. The optical path difference.
    :return: Float. The free spectral range.
    """
    return const.c / opd


def thermal_expansion(alpha, t, t_ref):
    """
    Calculate the expansion of a material, (1 + alpha * dt), based on L = L0(1 + alpha * dt).
    :param alpha: Float. The linear thermal expansion coefficient of the material.
    :param dt: Float. The change in temperature from the reference temperature (often 20 Celsius).
    :return: Float. The expansion of the material.
    """
    dt = t - t_ref
    return 1 + alpha * dt


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


def n_air(lam, t, p):
    """
    Calculate the index of refraction of air.

    :param lam: Float. The wavelength at which to calculate the refractive index, in microns.
    :param t: Float. The temperature, in celsius, at which to calculate the refractive index.
    :param p: Float. The pressure, in atmospheres, at which to calculate the refractive index.
    :return: Float. The refractive index of air.
    """
    return 1 + (6432.8 + 2949810 * lam**2 / (146 * lam**2 - 1) + 25540 * lam**2 / (41 * lam**2 - 1)) / (1 + (t - 15) * 3.4785e-3) * 1e-8 * p


def snell_cos_to_sin(n, theta):
    """
    Go from cos to sin using Snell's law.

    For the FWMI, one of the materials is always going to be air.
    :param n: Float. The index of refraction.
    :param theta: Float. The angle of incidence on the FWMI.
    :return: Float. Cos in terms of n and sin.
    """
    return np.sqrt(1 - np.sin(theta)**2 / n**2)


def opd_exact(n1, d1, n2, d2, n3, d3, theta):
    """
    Calculate the exact optical path difference for the hybrid FWMI.

    Eq. 2 in [1].

    When theta = theta_t, the tilt angle, then this function gives the fixed optical path difference (FOPD).
    :param n1: Float. The index of refraction for the pure glass arm.
    :param d1: Float. The thickness of the pure glass arm.
    :param n2: Float. The index of refraction for the glass part of the hybrid arm.
    :param d2: Float. The thickness of the glass part of the hybrid arm.
    :param n3: Float. The index of refraction for the air part of the hybrid arm.
    :param d3: Float. The thickness of the air part of the hybrid arm.
    :param theta: Float. The incident angle on the FWMI.
    :return: Float. The exact optical path difference.
    """
    return 2 * (n1 * d1 * snell_cos_to_sin(n1, theta)
                - n2 * d2 * snell_cos_to_sin(n2, theta)
                - n3 * d3 * snell_cos_to_sin(n3, theta))


def transmittance_simple(gamma, fopd):
    """
    Calculate the FWMI transmittance with no tilt.

    Eq. 24 in [1].
    :param gamma: Float. The spectral width of the incident signal.
    :param fopd: Float. The fixed optical path difference.
    :return: Float. The transmittance of the FWMI.
    """
    return 1 / 2 - 1 / 2 * np.exp(-(np.pi * gamma / (const.c / fopd))**2)


def sdr_simple(gamma_m, gamma_a, fopd):
    """
    Calculate the spectral discrimination ratio with no tilt.

    A version of Eq. 17 in [1].
    :param gamma_m: Float. The spectral width of the molecular signal.
    :param gamma_a: Float. The spectral width of the aerosol signal.
    :param fopd: Float. The fixed optical path difference.
    :return: Float. The spectral discrimination ratio.
    """
    return transmittance_simple(gamma_m, fopd) / transmittance_simple(gamma_a, fopd)


def delta_phi(opd_theta, opd_theta_t, nu):
    """
    Calculate the phase difference between a ray incident at angle theta and a ray at angle theta_t.

    :param opd_theta: Float. The optical path difference for the non-central incident angle.
    :param opd_theta_t: Float. The optical path difference for the central (tilt) incident angle.
    :param nu: Float. Center frequency of the transmitted laser.
    :return: Float. The phase difference.
    """
    return 2 * np.pi * nu * (opd_theta - opd_theta_t) / const.c


def local_transmittance(d_phi, i1, i2, gamma, fsr):
    """
    Calculate the local transmittance.

    :param d_phi: Float. The phase difference between a ray incident at some angle and a ray incident at the tilt angle.
    :param fsr: Float. The FSR of the FWMI. A function of the tilt angle.
    :param i1: Float. Intensity of the beam through arm 1.
    :param i2: Float. Intensity of the beam through arm 2.
    :param gamma: Float. The linewidth of the signal (molecular or aerosol).
    :param fsr: Float. The FSR of the FWMI. A function of the tilt angle.
    :return: Float. The local transmittance.
    """
    return i1 + i2 - 2 * np.sqrt(i1 * i2) * np.exp(-(np.pi * gamma / fsr)**2) * np.cos(d_phi)


def mapping_angle(rho, phi, theta_t, f):
    """
    Calculate the angle used in the transmittance map function.

    Eq. 16 in [1].
    :param rho: Float. The radial coordinate in the exit lens's image plane.
    :param phi: Float. The azimuthal coordinate in the exit lens's image plane.
    :param theta_t: Float. The tilt angle.
    :param f: Float. The exit lens's focal length.
    :return: Float. The angle to use with the local transmittance mapping function.
    """
    return np.arccos((2 * f * np.cos(theta_t)**2 - rho * np.sin(2 * theta_t) * np.cos(phi)) / (2 * np.sqrt(f**2 + rho**2) * np.cos(theta_t)))


def transmittance_map(del_phi, li_args):
    """
    Calculate the map from the local transmittance to the exit lens image plane.

    :param del_phi: Float. The phase difference between a ray incident at angle theta and a ray at the tilt angle.
    :param li_args: List. The other arguments for the local transmittance function: i1, i2, gamma, fsr
    :return: Float. The transmittance map.
    """
    return local_transmittance(del_phi, *li_args)


def transmittance_map_integrand(rho, phi, del_phi, li_args):
    """
    Calculate the map from the local transmittance to the exit lens image plane.

    :param rho: Float. The radial coordinate in the exit lens's image plane.
    :param del_phi: Float. The phase difference between a ray incident at angle theta and a ray at the tilt angle.
    :param li_args: List. The other arguments for the local transmittance function: i1, i2, gamma, fsr
    :return: Float. The integrand of the overall transmittance.
    """
    return rho * local_transmittance(del_phi, *li_args)


def overall_transmittance(l_map_args, theta_d, f):
    """

    :param l_map_args: List. The arguments for transmittance_map_integrand: theta_t, f, nu, li_args
    :param theta_d: Float. The half divergent angle that sets the limit on the integrand.
    :return:
    """
    integral = nquad(transmittance_map_integrand, [[0, f * theta_d], [-np.pi, np.pi]], l_map_args)[0]
    return integral / (np.pi * f**2 * theta_d**2)
