import math
import numpy as np
from viktor.errors import UserError, InputViolation
from valueengineering.material import *


def get_input_error(params, **kwargs):
    input_fields = dict(P=params.P, L=params.L, D=params.D, B=params.B, gc=params.gc, gs=params.gs,
                        h=params.h, c=params.c, a=params.a, lb_rqd=params.lb_rqd)
    violations = []
    for key, value in input_fields.items():
        if value is None:
            violations.append(InputViolation("Type a value", fields=[key]))

    # ds = float(params.ds_str[1:])
    # if params.c > params.h/2 - ds:
    #     violations.append(InputViolation("h >= 2*(c+ø_s) ikke overholdt", fields=["h", "c", "ds_str"]))

    if violations:
        raise UserError("You must fix invalid field values before a result can be obtained.",
                        input_violations=violations)


def get_R1(params, **kwargs):
    L = params.L / 1000  # fundamentsbredde, m
    R1 = (L**2 / math.pi)**.5

    return R1


def out_q_round(params, **kwargs):
    try:
        P = params.P            # Søjlelast, kN
        L = params.L / 1000     # fundamentsbredde, m
        q = P / L**2
        return round(q)
    except:
        return "N/A"


def get_R0(params, **kwargs):
    if params.DB_str == 'Cirkulær':
        D = params.D / 1000  # Diameter af søjle, m
        R0 = D / 2  # Radius af søjle, m
    else:
        B = params.B / 1000  # Bredde af søjle, m
        R0 = (B**2 / math.pi)**.5  # Radius af søjle, m

    return R0


def get_R(params, **kwargs):
    R0 = get_R0(params)  # Radius af søjle, m
    R1 = get_R1(params)  # Radius af fundament, m
    R = (R0 * R1 ** 2) ** (1 / 3)  # Radius ved overgang konstant til aftagende moment m_t

    return R


def get_m(params, **kwargs):
    P = params.P  # Søjlelast, kN

    R1 = get_R1(params)  # Radius af fundament, m
    R = get_R(params)  # Radius ved overgang konstant til aftagende moment m_t
    R0 = get_R0(params)  # Radius af søjle, m

    D = 2 * R0  # Diameter af søjle (Equivalent diameter ved kvadratisk søjle), m
    n_step = 100

    r_lst1 = list(np.linspace(0, R0, num=n_step, endpoint=False))
    r_lst2 = list(np.linspace(R0, R, num=n_step, endpoint=False))
    r_lst3 = list(np.linspace(R, R1, num=n_step))

    r_lst = r_lst1 + r_lst2 + r_lst3

    mt_lst = []
    mr_lst = []

    for r in r_lst1:
        mt = P / (2 * math.pi) * (1 - (D ** 2 / (4 * R1 ** 2)) ** (1 / 3))
        mr = P / (2 * math.pi) * (
                1 / 3 * (1 / (R1 ** 2) - 1 / ((0.5 * D) ** 2)) * r ** 2 + 1 - (D ** 2 / (4 * R1 ** 2)) ** (1 / 3))

        mt_lst.append(mt)
        mr_lst.append(mr)

    for r in r_lst2:
        mt = P / (2 * math.pi) * (1 - (D ** 2 / (4 * R1 ** 2)) ** (1 / 3))
        mr = P / (2 * math.pi) * \
             (1 / 3 * (r / R1) ** 2 - (D ** 2 / (4 * R1 ** 2)) ** (1 / 3) - (
                     1 / 3 * (R / R1) ** 2 - (D ** 2 / (4 * R1 ** 2)) ** (1 / 3)) * (R / r))

        mt_lst.append(mt)
        mr_lst.append(mr)

    for r in r_lst3:
        mt = P / (2 * math.pi) * (1 - (r / R1) ** 2)
        mr = 0

        mt_lst.append(mt)
        mr_lst.append(mr)

    return [r_lst, mt_lst, mr_lst]


def out_m_max(params, **kwargs):

    try:
        [r_lst, mt_lst, mr_lst] = get_m(params)
        m_max = max(mt_lst + mr_lst)
        m_max_round = round(m_max, 1)
    except:
        m_max_round = "N/A"

    return m_max_round


def get_mR(params, **kwargs):
    h = params.h  # mm
    c = params.c  # mm
    ds = float(params.ds_str[1:])  # mm
    a = params.a  # mm

    lb_rqd = params.lb_rqd / 1000  # Forankring af armering uden hensyntagen til opbuk, m
    R1 = get_R1(params)  # Radius af fundament, m

    # Materialeparametre
    fcd = get_concrete_params(params.fc_str, params.gc)['fcd']  # MPa
    fbd = get_concrete_params(params.fc_str, params.gc)['fbd']  # MPa
    fyk = get_rebar_params(params.YK_str, params.gs)['fyk']  # MPa
    fyd = get_rebar_params(params.YK_str, params.gs)['fyd']  # MPa

    # Fuld forankringslængde
    lb = fyk * ds / (4 * fbd)  # mm (forudsætter gode forankringsforhold (bunden af et fundament)).
    lb = lb / 1000  # m
    rb = max(min(R1 + lb_rqd - lb, R1), 0)  # radius ud til fuld forankring start, m
    r2_lst = [-R1, -R1, -rb, rb, R1, R1]

    # Momentbæreevne plade
    mR = get_mRd_rc_plate(h, c, fcd, ds, a, fyd)  # kNm/m
    mR_edge = mR * lb_rqd / lb
    mR_max = min(mR, mR * (R1 + lb_rqd) / lb)
    mR_lst = [0, mR_edge, mR_max, mR_max, mR_edge, 0]

    return [[rb, mR_max], [R1, mR_edge]]


def out_OK(params, **kwargs):
    try:
        [r_lst, mt_lst, mr_lst] = get_m(params)
        [r2_lst, mR_lst] = get_mR(params)
    except:
        return ""

    r1 = r2_lst[-3]
    mR1 = mR_lst[-3]
    r2 = r2_lst[-2]
    mR2 = mR_lst[-2]

    if mt_lst[0] > mR1:
        return "Bæreevne ikke OK!"

    for i, r in enumerate(r_lst):
        if r < r1:
            continue
        else:
            mR = mR1 + (mR2 - mR1) / (r2 - r1) * (r - r1)

            if round(mt_lst[i]) > round(mR):
                return "Bæreevne ikke OK!"

    return "Bæreevne OK!"


def vis_D(params, **kwargs):
    return params.DB_str == 'Cirkulær'


def vis_B(params, **kwargs):
    return params.DB_str == 'Kvadratisk'


def get_vol_c(params, **kwargs):
    L = params.L / 1000 # m
    h = params.h / 1000 # m
    vol_c = L**2 * h
    return vol_c


def out_vol_c(params, **kwargs):
    try:
        vol_c_round = round(get_vol_c(params), 2)
        return vol_c_round
    except:
        return "N/A"


def get_mass_s(params, **kwargs):
    rho_s = 7850            # kg/m3
    L = params.L / 1000     # m
    h = params.h / 1000     # m
    ds = float(params.ds_str[1:]) / 1000    # m
    a = params.a / 1000     # m
    As = math.pi / 4 * ds ** 2 / a          # m2/m

    vol_s = 2 * As * L
    mass_s = vol_s * rho_s
    return mass_s


def out_mass_s(params, **kwargs):
    try:
        mass_s_round = round(get_mass_s(params), 1)
        return mass_s_round
    except:
        return "N/A"