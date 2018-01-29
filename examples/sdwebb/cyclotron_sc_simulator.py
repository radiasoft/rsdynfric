import numpy as np
from matplotlib import pyplot as plt

# Particle initial conditions
v_x = 100.  #cm/sec
v_y = -57.  #cm/sec

x = .01  #cm
y = 0.   #cm

# Plasma densities
n_beam = 1.e15  #cm^-3
n_elec = 1.e15  #cm^-3

# B_field
B0 = 1.  #Tesla

# number of cyclotron orbits to simulate
n_cycles = 100

#
# Derived quantities for simulating things
#

e_charge = 4.80320427e-10     #esu
speed_of_light = 29979245800  #cm/sec
e_mass = 9.10938215e-28       #grams

B0 *= 10000. # Tesla to Gauss
omega_c = e_charge * B0 / (e_mass * speed_of_light)
omega_beam_sqrd = 4.*np.pi * n_beam * e_charge * e_charge / e_mass
omega_elec_sqrd = 4.*np.pi * n_elec * e_charge * e_charge / e_mass
delta_sqrd = omega_beam_sqrd - omega_elec_sqrd

# resolve the hell out of the cyclotron period

steps_per_period = 20
dt = (2.*np.pi/omega_c) / steps_per_period
n_steps = n_cycles*steps_per_period

# Compute the initial canonical momenta

px = e_mass * v_x - e_charge*B0*y/(speed_of_light)
py = e_mass * v_y

#
# Define some symplectic integrator functions
#

def radial_kick(x, y):
	return -delta_sqrd*x, -delta_sqrd*y

def drift(p):
	return p/e_mass

def similarity(x, y):
	return -e_mass*omega_c*y, -e_mass*omega_c*x

def similarity_inverse(x, y):
	return e_mass*omega_c*y, e_mass*omega_c*x


z_array = np.zeros([4,n_steps])

step = 0

while step < n_steps:

	z_array[0,step] = px
	z_array[1,step] = x
	z_array[2,step] = py
	z_array[3,step] = y

	px += - 0.5*dt * radial_kick(x, y)[0]
	py += - 0.5*dt * radial_kick(x, y)[1]
	
	y += 0.5*dt*drift(py)
	
	px += similarity(x, y)[0]
	py += similarity(x, y)[1]
	x += dt*drift(px)
	px += similarity_inverse(x, y)[0]
	py += similarity_inverse(x, y)[1]
	y += 0.5*dt*drift(py)

	px += - 0.5*dt * radial_kick(x, y)[0]
	py += - 0.5*dt * radial_kick(x, y)[1]

	step += 1

# Plot x versus y
plt.plot(z_array[1,:], z_array[3,:])
plt.show()


