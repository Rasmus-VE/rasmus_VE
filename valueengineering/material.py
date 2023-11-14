import numpy as np


def get_concrete_params(fc_str, gc):
    fck = float(fc_str[1:])
    fctm = 0.3 * fck ** (2 / 3)
    fctk = 0.7 * fctm
    fcd = fck / gc
    fctd = fctk / gc
    fbd = 2.25 * fctd

    return {'fck': fck, 'fctm': fctm, 'fctk': fctk, 'fcd': fcd, 'fctd': fctd, 'fbd': fbd}


def get_rebar_params(YK_str, gs):
    fyk = 500 if YK_str == 'K' else 550  # MPa
    fyd = fyk / gs

    return {'fyk': fyk, 'fyd': fyd}


def get_mRd_rc_plate(h, c, fcd, ds, a, fyd):
    As = np.pi / 4 * ds ** 2 / a  # mm2/mm
    d = h - (c + ds)  # mm
    omega = As * fyd / (d * fcd)  # enhedsløs
    mu = omega * (1 - omega / 2)  # enhedsløs
    mR = mu * d ** 2 * fcd  # momentbæreevne, Nmm/mm
    mR = mR / 1000  # momentbæreevne, kNm/m

    return mR
