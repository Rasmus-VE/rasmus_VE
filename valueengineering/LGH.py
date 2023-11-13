import math
import numpy as np
from viktor.errors import UserError, InputViolation


def get_input_error(params, **kwargs):
    input_fields = dict(P=params.P, L=params.L, D=params.D, B=params.B, gc=params.gc, gs=params.gs,
                        h=params.h, c=params.c, a=params.a, lb_rqd=params.lb_rqd)
    violations = []
    for key, value in input_fields.items():
        if value is None:
            violations.append(InputViolation("Type a value", fields=[key]))

    ds = float(params.ds_str[1:])
    if params.c > params.h/2 - ds:
        violations.append(InputViolation("h >= 2*(c+ø_s) ikke overholdt", fields=["h", "c", "ds_str"]))

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
    if params.DB == 'Cirkulær':
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

    r_lst1 = np.linspace(0, R0, num=n_step, endpoint=False)
    r_lst2 = np.linspace(R0, R, num=n_step, endpoint=False)
    r_lst3 = np.linspace(R, R1, num=n_step)

    r_lst = r_lst1.tolist() + r_lst2.tolist() + r_lst3.tolist()

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

    fck = float(params.fc_str[1:])
    fctm = 0.3 * fck ** (2 / 3)
    fctk = 0.7 * fctm
    gc = params.gc
    fcd = fck / gc
    fctd = fctk / gc
    fbd = 2.25 * fctd

    ds = float(params.ds_str[1:])  # mm
    a = params.a  # mm
    As = math.pi / 4 * ds ** 2 / a  # mm2/mm

    fyk = 500 if params.st_YK == 'K' else 550  # MPa
    gs = params.gs
    fyd = fyk / gs

    d = h - (c + ds)  # mm
    omega = As * fyd / (d * fcd)  # enhedsløs
    mu = omega * (1 - omega / 2)  # enhedsløs
    mR = mu * d ** 2 * fcd  # momentbæreevne, Nmm/mm
    mR = mR / 1000  # momentbæreevne, kNm/m

    lb = fyk * ds / (4 * fbd)  # mm (forudsætter gode forankringsforhold (bunden af et fundament)).
    lb = lb / 1000  # Fuld forankringslængde m
    lb_rqd = params.lb_rqd / 1000  # Forankring af armering uden hensyntagen til opbuk, m
    R1 = get_R1(params)  # Radius af fundament, m

    if lb < lb_rqd:
        mR_lst = [0, mR, mR, mR, mR, 0]
        r_lst2 = [-R1, -R1, 0, R1, R1]
    elif lb < R1:
        mR_lst = [0, mR * lb_rqd / lb, mR, mR, mR * lb_rqd / lb, 0]
        r_lst2 = [-R1, -R1, -(R1 - (lb - lb_rqd)), R1 - (lb - lb_rqd), R1, R1]
    else:
        mR_lst = [0, mR * lb_rqd / lb, mR * R1 / lb, mR * lb_rqd / lb, 0]
        r_lst2 = [-R1, -R1, 0, R1, R1]

    return [mR_lst, r_lst2]


def out_OK(params, **kwargs):
    try:
        [r_lst, mt_lst, mr_lst] = get_m(params)
        [mR_lst, r_lst2] = get_mR(params)
    except:
        return ""

    r1 = r_lst2[-3]
    mR1 = mR_lst[-3]
    r2 = r_lst2[-2]
    mR2 = mR_lst[-2]
    # print("r1 =", r1, ", mR1 =", mR1, ", r2 =", r2, " mR2 =", mR2)

    if mt_lst[0] > mR1:
        return "Bæreevne ikke OK!"

    for i, r in enumerate(r_lst):
        if r < r1:
            continue
        else:
            mR = mR1 + (mR2 - mR1) / (r2 - r1) * (r - r1)
            # print("mR = ", mR, ", mE = ", mt_lst[i])

            if round(mt_lst[i]) > round(mR):
                return "Bæreevne ikke OK!"

    return "Bæreevne OK!"


def vis_D(params, **kwargs):
    return params.DB == 'Cirkulær'


def vis_B(params, **kwargs):
    return params.DB == 'Kvadratisk'


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
    print(As)
    vol_s = 2 * As * L
    mass_s = vol_s * rho_s
    return mass_s


def out_mass_s(params, **kwargs):
    try:
        mass_s_round = round(get_mass_s(params), 1)
        return mass_s_round
    except:
        return "N/A"