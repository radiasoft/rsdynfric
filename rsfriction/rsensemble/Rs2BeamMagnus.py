# -*- coding: utf-8 -*-

"""Encapsulation of a 3D electron bunch, with N ions.
We are considering the two-particle problem described by

D.L. Bruhwiler and S.D. Webb, “New algorithm for dynamical friction
of ions in a magnetized electron beam,” AIP Conf. Proc. 1812, 050006
(2017); https://doi.org/10.1063/1.4975867

The number of ions N will typically be 1 or a few. The number of
electrons must be large enough to obtain a converged result.

We are essentially doing a Monte Carlo integration of the two-body
problem. The ions do not see each other. Each ion sees the original
electron distribution -- not the perturbed electrons.

The ion and electron coordinates are specified by the user,
according to the Hamiltonian from Eqs. (1) and (2) of the paper.
We assume magnetized electrons, executing fast Larmor oscillations.
We assume drifting ions, for which the magnetic field is negligible.

However, all coordinates must be transformed into the adiabatically
averaged action-angle coordinates shown in Eqs. (9) and (10) of
the paper.

In these coordinates, the momentum change of each ion is calculated
for each electron. The sum of changes, divided by the interaction
time, yields the dynamic friction force that is exerted on each
ion. This is fundamental result to be calculated here.

For the 'original' phase space, we assume 3 phase space coordinates
and 3 corresponding momenta, with the following conventions:
    x [m], xp=px/p0 [rad]
    y [m], yp=py/p0 [rad]
    s [m], dp=(p-p0)/p0 [rad]

where s is the total distance traveled along the accelerator axis,
and p0 [eV/c] is the design momentum.

The physical problem being addressed here involves relativistic,
copropagating electron and ion beams. However, the dynamic friction
calculation takes place in the beam frame, where p0=0.

In the beam frame, all velocities are non-relativistic. In this frame,
the electron beam is a bounded, thermal plasma.

For the electron beam module, imported from the 'rsbeams' library,
  https://github.com/radiasoft/rsbeams/
it is required to specify all Twiss parameters. In this class,
however, we assume the beam is close to a waist, transversely and
longitudinally, so that all Twiss alpha parameters are zero.

We also allow the user to specify RMS values of position and
momentum for the electron beam, and a beam shape, rather than
having to specify emittances and Twiss beta values.

Ion phase space coordinates and momenta are specified explicitly
for each particle, and they must be normalized to the RMS value
of the corresponding coordinate/momentum for the e- distribution.

:copyright: Copyright (c) 2018 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import math
import numpy

from rsbeams.rsstats import stats6d
from rsbeams.rsptcls import RsTwiss2D
from rsbeams.rsptcls import RsPtclBeam6D
from rsbeams.rsphysics import rsconst

class Rs2BeamMagnus:
    """Representation of a 6D charged particle distribution."""

    def __init__(self, n_, design_p_ev, total_charge_c, mass_ev, dist_type, max_rms_fac, \
                 alpha_x, beta_x, emit_x, alpha_y, beta_y, emit_y, alpha_z, beta_z, emit_z):

        x = 10.


    def test_beam_gen():

        # specify physical properties of the beam
        num_ptcls = 1000
        design_p_ev = 271e+6
        total_charge_c = 3.05e-09
        mass_ev = rsconst.m_e_EV

        dist_type = 'gaussian'
        max_rms_fac = 4.9

        alpha_x = 1.3     # []
        beta_x = 21.1     # [m/rad]
        emit_x = 1.7e-06  # [m-rad]

        alpha_y = -3.01   # []
        beta_y = 9.07     # [m/rad]
        emit_y = 0.44e-6  # [m-rad]

        alpha_z = 0       # []
        beta_z =  1.55    # [m/rad]
        emit_z = 3.1e-05  # [m-rad]

        my_ebeam = RsPtclBeam6D.RsPtclBeam6D(num_ptcls, design_p_ev, \
                                             total_charge_c, mass_ev, \
                                             dist_type, max_rms_fac, \
                                             alpha_x, beta_x, emit_x, \
                                             alpha_y, beta_y, emit_y, \
                                             alpha_z, beta_z, emit_z )

        beta_gamma = my_ebeam.get_beta0_gamma0()
        gamma0 = my_ebeam.get_gamma0()
        beta0 = my_ebeam.get_beta0()

#        print('beta_gamma = ', beta_gamma)
#        print('gamma0     = ', gamma0)
#        print('beta0      = ', beta0)

        my_twiss_x = my_ebeam.get_twiss2d_by_name('twiss_x')
        my_twiss_y = my_ebeam.get_twiss2d_by_name('twiss_y')
        my_twiss_z = my_ebeam.get_twiss2d_by_name('twiss_z')

#        print('alpha_x = ', my_twiss_x.get_alpha_rms())
#        print('beta_x = ', my_twiss_x.get_beta_rms())
#        print('emit_x = ', my_twiss_x.get_emit_rms())
