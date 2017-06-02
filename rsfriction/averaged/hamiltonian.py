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

def _calc_omega_larmor(b_field, mass_number=1, charge_number=-1):
    """Calculate the angular frequency of Larmor rotations, which is a signed quantity.
    
    Args:
        b_field:        uniform magnetic field amplitude [T]
        mass_number:    charged particle mass, normalized to electron [kg]
        charge_number:  charge number of the particle []
        
    Returns:
        omega_larmor:   the Larmor frequency [radians/s]
    """

    omega_larmor = (charge_number * elementary_charge) * b_field / (mass_number * electron_mass)
    return omega_larmor

def _calc_hamiltonian_0(z_i, z_e, ion_mass, b_field):
    """Calculate the Hamiltonian for a drifting ion and a magnetized electron.
   
    Args:
        z_i: 6xN array of ion coördinates and conjugate momenta
        z_e: 6xN array of electron coördinates and conjugate momenta
             Phase-space variables are (x, px, y, py, z, pz)
        ion_mass:  mass of the ion [Kg]
        b_field:   uniform magnetic field amplitude [T]

    Returns:
        hamiltonian_0: the lowest-order Hamiltonian of each ion-electron pair
    """

    h_ion = ((z_i[i_px,:]**2 + z_i[i_py, :]**2 + z_i[i_pz, :]**2) / (2.*ion_mass))
    h_electron = ((z_e[i_px,:] - elementary_charge*b_field*z_e[i_y,:])**2 + z_e[i_py,:]**2 + z_e[i_pz,:]**2) / (2.*electron_mass)

    hamiltonian_0 = h_ion + h_electron
    return hamiltonian_0

def _calc_hamiltonian_c(z_i, z_e, ion_charge):
    """Calculate the Hamiltonian for the Coulomb interaction of each ion-electron pair.
    
    Args:
        z_i: 6xN array of ion coördinates and conjugate momenta
        z_e: 6xN array of electron coördinates and conjugate momenta
             Phase-space variables are (x, px, y, py, z, pz)
        ion_mass:  mass of the ion [Kg]
        b_field:   uniform magnetic field amplitude [T]

    Returns:
        hamiltonian_0: the lowest-order Hamiltonian of each ion-electron pair
    """
    g_ie = -ion_charge*elementary_charge / (4.*pi*epsilon_0)
    hamiltonian_c = g_ie / np.sqrt(
        + (z_i[i_x, :] - z_e[i_x, :])**2
        + (z_i[i_y, :] - z_e[i_y, :])**2
        + (z_i[i_z, :] - z_e[i_z, :])**2)
    return hamiltonian_c

def _calc_hamiltonian(z_i, z_e, ion_mass, ion_charge, b_field):
    """Calculate the total Hamiltonian for each ion-electron pair.
    """
    hamiltonian = _calc_hamiltonian_0(z_i, z_e, ion_mass, b_field)\
                  + _calc_hamiltonian_c(z_i, z_e, ion_charge)
    return hamiltonian
