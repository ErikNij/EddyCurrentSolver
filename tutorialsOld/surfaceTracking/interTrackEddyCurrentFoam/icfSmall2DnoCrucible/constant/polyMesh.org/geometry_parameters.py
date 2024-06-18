#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale = 1e-3

geo_xS =   0.0
geo_x0 =  85.0
geo_x1 = geo_x0 +  1.0
geo_x2 = geo_x1 + 14.0
geo_x3 = geo_x2 + 15.0
geo_x4 = 300.0

geo_z0 = geo_x1 - geo_x0
geo_z1 = 200.0
geo_z2 = 235.0
geo_z3 = 383.0
geo_z4 = geo_z1 + geo_z3

axi_xS = geo_xS
axi_x0 = geo_x0
axi_x1 = geo_x1
axi_x2 = geo_x2
axi_x3 = geo_x3
axi_x4 = geo_x4

n_scale = 0.25
#n_scale = 1.0

dy = 0.5 * 1.0/n_scale

axi_yS = dy
axi_y0 = dy
axi_y1 = dy
axi_y2 = dy
axi_y3 = dy
axi_y4 = dy

n_x0 = int(m.ceil(n_scale*geo_x0))
n_x1 = int(m.ceil(n_scale*(geo_x1-geo_x0)))
n_x2 = int(m.ceil(n_scale*(geo_x2-geo_x1)))
n_x3 = int(m.ceil(n_scale*(geo_x3-geo_x2)))
n_x4 = int(m.ceil(n_scale*(geo_x4-geo_x3)))

n_z0 = int(m.ceil(n_scale*geo_z0))
n_z1 = int(m.ceil(n_scale*(geo_z1-geo_z0)))
n_z2 = int(m.ceil(n_scale*geo_z2))
n_z3 = int(m.ceil(n_scale*(geo_z3-geo_z2)))
n_z4 = int(m.ceil(n_scale*(geo_z4-geo_z3)))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
