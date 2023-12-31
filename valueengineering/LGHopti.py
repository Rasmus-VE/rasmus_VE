import math
import numpy as np
import pandas as pd
from valueengineering.LGH import get_R1, get_R0, get_m, get_vol_c
from viktor.errors import InputViolation, UserError
from valueengineering.LGH import *
from valueengineering.material import *


def get_input_error_opti(params, **kwargs):
    input_fields = dict(P=params.P, L=params.L, D=params.D, B=params.B, gc=params.gc, gs=params.gs,
                        h=params.h, c=params.c, a_min=params.a_min, a_max=params.a_max, lb_rqd=params.lb_rqd)
    violations = []
    for key, value in input_fields.items():
        if value is None:
            violations.append(InputViolation("Type a value", fields=[key]))

    if violations:
        raise UserError("You must fix invalid field values before a result can be obtained.",
                        input_violations=violations)


def get_mR_opti(params, ds_str, a, **kwargs):
    """
    get_mR_opti returnerer koordinatsæt til plot af momentbæreevnen mR.
    """
    params.ds_str = ds_str
    params.a = a

    [[rb, mR_max], [R1, mR_edge]] = get_mR(params)

    r2_lst = [-R1, -R1, -rb, rb, R1, R1]
    mR_lst = [0, mR_edge, mR_max, mR_max, mR_edge, 0]

    return [r2_lst, mR_lst]


def get_ds_a_opti(params, **kwargs):
    """
    Sorterer dataframe med alle brugbare løsninger efter mass og udvælger den letteste (første i liste)
    """
    case_lst = get_case_lst(params)
    if not case_lst:
        violation = InputViolation("", fields=["a_min", "a_max", "ds_str_lst"])
        raise UserError("Ingen armering inden for de valgte grænser kan overholde bæreevnen",
                        input_violations=[violation])

    df = pd.DataFrame(data=case_lst)
    df = df.sort_values(by=['As'])

    return df.iloc[0, 0], df.iloc[0, 1]


def get_case_lst(params, **kwargs):
    """
    Beregner momentbæreevnen af alle kombonationer af ds og a. Frasorterer alle design der ikke for tilstrækkelig styrke.
    Returnerer alle OK kombinationer af ds, a og As.
    """
    [r_lst, mt_lst, mr_lst] = get_m(params)
    mE = max(mt_lst)

    ds_array = params.ds_str_lst  # mm

    a_min, a_max, a_step = params.a_min, params.a_max, params.a_step
    a_start = a_min if a_min % a_step == 0 else a_min + a_step - (a_min % a_step)
    a_stop = a_max + a_step if a_max % a_step == 0 else a_max
    a_array = range(a_start, a_stop, a_step)  # mm

    case_lst = []

    for ds_str in ds_array:
        for a in a_array:
            params.ds_str = ds_str
            params.a = a

            [[rb, mR_max], [R1, mR_edge]] = get_mR(params)

            if mR_max < mE:
                continue

            # Tjek:
            OK = True
            for r, mt in zip(r_lst, mt_lst):
                if r < rb:
                    continue
                else:
                    mR = mR_max + (mR_edge - mR_max) / (R1 - rb) * (r - rb)
                    if round(mt) > round(mR):
                        OK = False
            if OK is False:
                continue

            # Calculate As
            ds = float(params.ds_str[1:])
            As = math.pi / 4 * ds ** 2 / a  # mm2/mm

            case = {'ds': ds, 'a': a, 'As': As}
            case_lst.append(case)

    return case_lst


def get_mass_s_opti(params, **kwargs):
    rho_s = 7850            # kg/m3
    L = params.L / 1000     # m
    h = params.h / 1000     # m
    ds, a = get_ds_a_opti(params)   # mm, mm
    ds = ds / 1000    # m
    a = a / 1000     # m
    As = math.pi / 4 * ds ** 2 / a          # m2/m
    vol_s = 2 * As * L
    mass_s = vol_s * rho_s
    return mass_s


def out_mass_s_opti(params, **kwargs):
    try:
        mass_s_round = round(get_mass_s_opti(params), 1)
        return mass_s_round
    except:
        return "N/A"


def out_mass_s_per_m3_opti(params, **kwargs):
    try:
        mass_s_per_m3_opti = round(get_mass_s_opti(params) / get_vol_c(params), 1)
        return mass_s_per_m3_opti
    except:
        return "N/A"


def out_rebar_name(params, **kwargs):
    try:
        YK_str = params.YK_str
        ds, a = get_ds_a_opti(params)  # mm, mm
        mass_s = get_mass_s_opti(params)   # kg
        rebar_name = f" {YK_str}{str(int(ds))}/{a} us-br"
        return rebar_name
    except:
        return ""