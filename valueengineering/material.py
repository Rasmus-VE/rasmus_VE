
def get_concrete_params(fc_str, gc):
    fck = float(fc_str[1:])
    fctm = 0.3 * fck ** (2 / 3)
    fctk = 0.7 * fctm
    fcd = fck / gc
    fctd = fctk / gc
    fbd = 2.25 * fctd

    return {'fck': fck, 'fctm': fctm, 'fctk': fctk, 'fcd': fcd, 'fctd': fctd, 'fbd': fbd}
