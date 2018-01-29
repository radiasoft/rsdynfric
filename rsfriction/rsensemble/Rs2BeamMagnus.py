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

Unfortunately, p0=0 creates a problem for the particle beam classes
in the 'rsbeams' library, so we specify p0=1. and accommodate.

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
    """6D electron distribution plus a handful of arbitrary ions."""

    def __init__(self, num_elec, charge_ebeam_c, distrib_type, max_rms_fac, \
                 x_rms_m, vx_rms_mks, y_rms_m, vy_rms_mks, z_rms_m, vz_rms_mks, \
                 charge_ions_c, mass_ions_c):
        """Initialize electron beam and ion properties"""

        # Define the electron beam (p0=1 with Twiss alpha's set to zero)
        emit_x = math.sqrt(x_rms_m * vx_rms_mks) # [m-rad]
        beta_x = x_rms_m * (x_rms_m/emit_x)      # [m/rad]
        emit_y = math.sqrt(y_rms_m * vy_rms_mks) # [m-rad]
        beta_y = y_rms_m * (y_rms_m/emit_y)      # [m/rad]
        emit_z = math.sqrt(z_rms_m * vz_rms_mks) # [m-rad]
        beta_z = z_rms_m * (z_rms_m/emit_z)      # [m/rad]

        self.num_elec = 1000
        self.design_p_ev = 0.
        self.max_rms_fac = max_rms_fac
        self.mass_elec_ev = rsconst.m_e_EV
        self.total_charge_c = charge_ebeam_c

        if ( (distrib_type != 'uniform') and
             (distrib_type != 'gaussian') ):
            message = '\n\nERROR --'
            message += '\n    distrib_type is specified as "' + self.distrib_type + '", which is not supported.'
            message += '\n    Only "uniform" and "gaussian" are allowed.'
            message += '\n'
            raise Exception(message)
        else:
            self.distrib_type = distrib_type

        self.e_beam = RsPtclBeam6D.RsPtclBeam6D(self.num_elec, 0., self.total_charge_c, \
                                             self.mass_elec_ev, self.distrib_type, self.max_rms_fac, \
                                             0., beta_x, emit_x, \
                                             0., beta_y, emit_y, \
                                             0., beta_z, emit_z )

        my_twiss_x = self.e_beam.get_twiss2d_by_name('twiss_x')
        my_twiss_y = self.e_beam.get_twiss2d_by_name('twiss_y')
        my_twiss_z = self.e_beam.get_twiss2d_by_name('twiss_z')

#        print('alpha_x = ', my_twiss_x.get_alpha_rms())
#        print('beta_x = ', my_twiss_x.get_beta_rms())
#        print('emit_x = ', my_twiss_x.get_emit_rms())
