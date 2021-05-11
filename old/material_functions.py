from numpy import sqrt


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


def delta_n(t, n_ref, d0, d1, d2, e0, e1, lam_tk, lam, t_ref=20):
    """
    Calculate the change in refractive index, with respect to the reference index, as a function of temperature
    deviation from the reference temperature.

    This equation can be found in [1] as well as chapter 22 of the Zemax user manual.
    :param t: Float. The temperature at which to calculate the refractive index change.
    :param n_ref: Float. The reference refractive index at the reference temperature at the wavelength of interest.
    :param d0: Float.
    :param d1: Float.
    :param d2: Float.
    :param e0: Float.
    :param e1: Float.
    :param lam_tk: Float. Reference wavelength, in microns. Determined at the reference temperature, but is a fit
                   parameter and not an explicit function of index or temperature.
    :param lam: Float. Wavelength of interest, in microns.
    :param t_ref. Float. The temperature at which the reference refractive index was measured. Default is 20 Celsius.
    :return: Float. The temperature derivative of the refractive index.
    """
    dt = t - t_ref
    return (n_ref**2 - 1) / (2 * n_ref) * (d0 * dt + d1 * dt**2 + d2 * dt**3 + (e0 * dt + e1 * dt**2) / (lam**2 - lam_tk**2))


def beta(t, n_ref, d0, d1, d2, e0, e1, lam_tk, lam, t_ref=20):
    """
    Calculate the temperature derivative of the refractive index of glass.

    :param t: Float. The temperature at which to calculate the refractive index change.
    :param n_ref: Float. The reference refractive index at the reference temperature at the wavelength of interest.
    :param d0: Float.
    :param d1: Float.
    :param d2: Float.
    :param e0: Float.
    :param e1: Float.
    :param lam_tk: Float. Reference wavelength, in microns. Determined at the reference temperature, but is a fit
                   parameter and not an explicit function of index or temperature.
    :param lam: Float. Wavelength of interest, in microns.
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
