from __future__ import absolute_import, division, print_function, unicode_literals

import rsmath.const as rsc
from rsfriction.rsensemble import Rs2BeamMagnus
import scipy

def test_twobeam_magnus():

    # specify physical properties of the beam
    mass_ev = scipy.constants.m_e * rsc.C_SQ / scipy.constants.e  # [eV] mass of an electron
    num_elec = 1000
    max_rms_fac = 4.9
    distrib_type = 'gaussian'
    charge_ebunch_c = 1.0e-09  # [C] total charge of the electron distribution

    x_rms_m = 1.e-3    # [m]
    y_rms_m = 1.e-3    # [m]
    z_rms_m = 1.e-2    # [m]

    vx_rms_mks = 1.e+5    # [m/s]
    vy_rms_mks = 1.e+5    # [m/s]
    vz_rms_mks = 1.e+4    # [m/s]

    # specify physical properties of the ion(s)
    charge_ion_c = scipy.constants.e
    mass_ion_kg = scipy.constants.m_p   # [kg] mass of a proton
    mass_ion_ev = mass_ion_kg * rsc.KG_to_EV

    my_2beam = Rs2BeamMagnus.Rs2BeamMagnus(num_elec, charge_ebunch_c, distrib_type, max_rms_fac, \
                                           x_rms_m, vx_rms_mks, y_rms_m, vy_rms_mks, z_rms_m, vz_rms_mks, \
                                           charge_ion_c, mass_ion_ev)

#    print('alpha_y = ', new_twiss_y.get_alpha_rms())
#    print('beta_y = ', new_twiss_y.get_beta_rms())
#    print('emit_y = ', new_twiss_y.get_emit_rms())

#    assert(stats6d.specify_significant_figures(new_twiss_y.get_alpha_rms(),3) == \
#           stats6d.specify_significant_figures(alpha_y,3))
#    assert(stats6d.specify_significant_figures(new_twiss_y.get_beta_rms(),3) == \
#           stats6d.specify_significant_figures(beta_y,3))
#    assert(stats6d.specify_significant_figures(new_twiss_y.get_emit_rms(),3) == \
#           stats6d.specify_significant_figures(emit_y,3))

test_twobeam_magnus()
