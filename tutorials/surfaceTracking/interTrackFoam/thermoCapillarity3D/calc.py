#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

import math as m

wedge = 5.0

L = 4.00 * 1e-3               # [L]         = m
R = 0.15 * 1e-3               # [R]         = m
nu = 1e-6                     # [nu]        = m^2/s
rho = 998.00                  # [rho]       = kg/m^3
eta = nu*rho                  # [eta]       = m^2/s * kg/m^3 = kg/s/m
ddxT = 25.4e+3/500            # [ddxT]      = K/m
#ddxT = 3.23e+3            # [ddxT]      = K/m
#ddxT = 25.4e+4            # [ddxT]      = K/m
ddTSigma = 0.0274e-3          # [ddTSigma]  = N/m/K
alpha = 1.08e-9*2000          # [alpha]     = m^2/s

V_i = 4.0/3 * m.pi * R**3     # [V_i]       = m^3
V_i_wedge = wedge/360 * V_i   # [V_i_wedge] = m^3
rho_i = 1.161                 # [rho_i]     = kg/m^3
m_i = V_i * rho_i             # [m_i]       = kg
m_i_wedge = V_i_wedge * rho_i # [m_i_wedge] = kg

#

#R = 0.1e-2
#eta = 0.24
#rho = 0.5
#nu = eta/rho
#ddxT = 2.0/1e-2
#ddTSigma = 2
#cp = 1.0
#alpha = 0.24/(rho*cp)

#

u_0 = R*ddxT/eta * ddTSigma # [u_0]      = m*K/m * s*m/kg * kg*m/m/K/s^2 = m/s
u_inf = 0.5 * u_0

Re = u_inf*R / nu
Pr = nu/alpha
Mg = ddTSigma * R*(R*ddxT) / (eta*alpha)
#Mg = Re*Pr

print("R            = %g m = %g mm" % (R,R*1e+3))
print("nu           = %g m^2/s" % (nu))
print("rho          = %g kg/m^3" % (rho))
print("eta          = %g kg/s/m" % (eta))
print("ddxT         = %g K/m" % (ddxT))
print("ddTSigma     = %g N/m/K" % (ddTSigma))
print("alpha        = %g m^2/s" % (alpha))


print("u_0          = %g m/s = %g mm/s" % (u_0,u_0*1e+3))
print("u_inf        = %g m/s = %g mm/s" % (u_inf,u_inf*1e+3))

print("Re           = %g" % (Re))
print("Pr           = %g" % (Pr))
print("Mg           = %g" % (Mg))

print("V_i          = %g m^3 = %g mm^3" % (V_i,V_i*1e+9))
print("V_i_wedge    = %g m^3 = %g mm^3" % (V_i_wedge,V_i_wedge*1e+9))
print("rho_i        = %g kg/m^3" % (rho_i))
print("m_i          = %g kg = %g g" % (m_i,m_i*1e+3))
print("m_i_wedge    = %g kg = %g g" % (m_i_wedge,m_i_wedge*1e+3))