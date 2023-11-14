from LGH import *
from munch import Munch

params = Munch()

params.DB = 'Rektangulær'
params.D = 500
params.B = 400
params.h = 450
params.c = 50
params.a = 100
params.ds_str = 'ø12'
params.gc = 1.45
params.fc_str = 'C30'
params.YK_str = 'Y'
params.P = 1200  # kN
params.L = 2800  # mm
params.lb_rqd = 0
params.YK_str = 'Y'
params.gs = 1.2
params.gc = 1.45

[r, mR] = get_mR(params)

print(r)
print(mR)
