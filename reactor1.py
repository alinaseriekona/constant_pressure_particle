"""
Constant-pressure, adiabatic kinetics simulation.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import sys

import cantera as ct
from Particle import *
ct.add_directory('./data')
gas = ct.Solution('caltech.yaml')
gas.TPX = 1700.0, 12*ct.one_atm, 'CH4:1,N2:1'
r = ct.IdealGasConstPressureReactor(gas, energy = 'off')

sim = ct.ReactorNet([r])
sim.verbose = True

# Particle arrays initilization
PND = [] # Particle Number Density
PND.append(0.0)
FV = [] # Particle Volume Fraction
FV.append(0.0)
X_A4R5C = [] # Mole fraction of A4R5 after consumption with the particle module
X_A4R5C.append(0.0)


# limit advance when temperature difference is exceeded
delta_T_max = 20.
r.set_advance_limit('temperature', delta_T_max)

dt_max = 1.e-4
t_end = 1000 * dt_max
states = ct.SolutionArray(gas, extra=['t'])

print('{:10s} {:10s} {:10s} {:14s}'.format(
    't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    states.append(r.thermo.state, t=sim.time*1e3)

    #******************************************************************************************
    # Particle formation calculations
    # 1. Calculate the inception rate and production of new particles
    particle_total_number = PND[-1] + inception_rate(X_A4R5C[-1],r.thermo.P,r.thermo.T)*36/500*dt_max 
    PND.append(particle_total_number)
    # 2. Calculate the volume fraction of the particles
    FV.append(volume_fraction(PND[-1]))
    # 3. Estimate the mole fraction of A4R5 after consumption by particle module
    X_A4R5C.append(gas_scrub(PND[-1],r.thermo.X[gas.species_index('A4R5')], r.thermo.P, r.thermo.T)) 


    #???????????
    # How to overwrite A4R5C on A4R5 before going to next iteration? How to update the composition using A4R5C?
    # If I can apply this change, then in the 1st step of particle formation, instead of using X_A4R5C[-1], I can
    # use r.thermo.X[gas.species_index('A4R5')]. The benefit of this mothod is that A4R5 consumption will be 
    # reflected on chemistry and move the equilibrium forward to form more A4R5.  
    #******************************************************************************************
    print('{:10.3e} {:10.3f} {:10.3f} {:14.6f} {:10.3e}'.format(sim.time, r.T, r.thermo.P, r.thermo.u, PND[-1]))

# Plot the results if matplotlib is installed.
# See http://matplotlib.org/ to get it.

if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(states.t, FV[:-1])
    plt.xlabel('Time (ms)')
    plt.ylabel('FV (ppm)')
    plt.subplot(2, 2, 2)
    plt.plot(states.t, states.X[:, gas.species_index('CH4')])
    plt.xlabel('Time (ms)')
    plt.ylabel('CH4 Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(states.t, X_A4R5C[:-1])
    plt.xlabel('Time (ms)')
    plt.ylabel('A4R5C Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(states.t, states.X[:, gas.species_index('A4R5')])
    plt.xlabel('Time (ms)')
    plt.ylabel('A4R5 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")
