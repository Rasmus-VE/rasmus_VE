import math

def get_R1(params, **kwargs):
    L = params.L / 1000  # fundamentsbredde, m
    R1 = (L**2 / math.pi)**.5

    return R1