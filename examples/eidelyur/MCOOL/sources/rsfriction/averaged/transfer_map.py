# -*- coding: utf-8 -*-
u"""Library of simple methods for charged particle dynamics.
:copyright: Copyright (c) 2017 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function
from scipy.constants import elementary_charge
from scipy.constants import electron_mass
from scipy.constants import epsilon_0
from scipy.constants import pi

# Specify indexing of the particle phase space coordinates
(i_x, i_px, i_y, i_py, i_z, i_pz) = range(6)

"""
define transfer maps for ions and electrons
There are three maps to define here: one each
for ions and electrons under H_0, and another
"""

""" matrix for a linear drift """
def MatD(mass, h):
    Mdrift = np.identity(6)
    for i in (Ix, Iy, Iz):
        Mdrift[i, i + 1] = h / mass
    return Mdrift

""" matrix for linear electron dynamics in a solenoidal field """
def MatK0_e(h):
    mw = me * Omega_e
    wh = Omega_e * h
    cwh = math.cos(wh)
    swh = math.sin(wh)
    cwh1m = 2 * math.sin(wh / 2) ** 2  # 1 - cos(a) = 2 sin^2(a / 2)
    MK0 = np.identity(6)
    MK0[Iy,  Iy ] = cwh
    MK0[Ipy, Ipy] = cwh
    MK0[Iy,  Ipy] = swh / mw
    MK0[Ipy, Iy ] = -mw * swh
    MK0[Iz,  Ipz] = h / me
    MK0[Ix,  Ipx] = swh / mw
    MK0[Ix,  Iy ] = swh
    MK0[Ix,  Ipy] = cwh1m / mw
    MK0[Iy,  Ipx] = -cwh1m / mw
    MK0[Ipy, Ipx] = -swh
    return MK0

"""
map phase-space coördinates forward in time by amount h
based on the Hamiltonian H_0, which describes the free
motion of ions and the motion of electrons in a solenoidal
magnetic field
"""
def MapZ_0(h, z_i, z_e):
    mat = MatD(M_ion, h)
    zf_i = mat.dot(z_i)
    mat = MatK0_e(h)
    zf_e = mat.dot(z_e)
    return zf_i, zf_e

"""
map phase-space coördinates forward in time by amount h
based on the Hamiltonian H_C, which describes the collision
between a single ion-electron pair
"""
def MapZ_C(h, z_i, z_e):
    g = h * Z_ion * qe ** 2 / (4 * pi * eps0)
    dz = z_i - z_e
    denom = (dz[Ix,:] ** 2 + dz[Iy,:] ** 2 + dz[Iz,:] ** 2) ** (3/2)
    zf_i = z_i.copy()
    zf_e = z_e.copy()
    for ip in (Ipx, Ipy, Ipz):
        zf_i[ip,:] = z_i[ip,:] - g * dz[ip - 1] / denom
        zf_e[ip,:] = z_e[ip,:] + g * dz[ip - 1] / denom
    return zf_i, zf_e

def apply_MapZ_0(h, n, z_i, z_e):
    mat_i = MatD(M_ion, h)
    mat_e = MatK0_e(h)
    zf_i = [z_i]
    zf_e = [z_e]
    for i in range(n):
        z_i = mat_i.dot(z_i)
        z_e = mat_e.dot(z_e)
        zf_i.append(z_i)
        zf_e.append(z_e)
    return np.asarray(zf_i), np.asarray(zf_e)

""" second-order split-operator integration for the total Hamiltonian """
def apply_MapZ(h, n, z_i, z_e):
    hh = 0.5 * h
    mat_i = MatD(M_ion, hh)
    mat_e = MatK0_e(hh)
    zf_i = [z_i]
    zf_e = [z_e]
    for i in range(n):
        z_i = mat_i.dot(z_i)
        z_e = mat_e.dot(z_e)
        z_i, z_e = MapZ_C(h, z_i, z_e)
        z_e = mat_e.dot(z_e)
        z_i = mat_i.dot(z_i)
        zf_i.append(z_i)
        zf_e.append(z_e)
return np.asarray(zf_i), np.asarray(zf_e)