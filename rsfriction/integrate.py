# -*- coding: utf-8 -*-
u"""Integrate an ion and a magnetized electron forward in time,
    using the original phase space coordinates.
    The fast Larmor oscillations of the electron must be resolved.

:copyright: Copyright (c) 2017 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function
import numpy
import math

from scipy.constants import pi
from scipy.constants import epsilon_0
from scipy.constants import proton_mass
from scipy.constants import electron_mass
from scipy.constants import Boltzmann
from scipy.constants import elementary_charge

# from scipy.constants import speed_of_light
# from scipy.constants import mu_0

from rsfriction.original import hamiltonian

fourPiEps0 = 4. * pi * epsilon_0
invFourPiEps0 = 1. / fourPiEps0

""" reset some default options """
# numpy.set_printoptions(linewidth=96)

""" prefixes """
(femto, pico, nano, micro, milli, one, kilo, mega, giga, tera, peta) = \
    10. ** numpy.asarray(range(-15, 15+1, 3))

z_ion = 1
m_ion = proton_mass
magnetic_field = 1. # [T] along z-axis
electron_temp = 300. # [deg K]
electron_charge = -1. * elementary_charge

n_gyro = 100 # a somewhat arbitrary choice, range [100, 160]

""" derived quantities """
v_th = math.sqrt(2. * Boltzmann * electron_temp / electron_mass)
rho_gc = electron_mass * v_th / (elementary_charge * magnetic_field)
omega_e = hamiltonian._calc_omega_larmor(magnetic_field)
T_e = (2 * pi) / abs(omega_e)
T_intxn = n_gyro * T_e

