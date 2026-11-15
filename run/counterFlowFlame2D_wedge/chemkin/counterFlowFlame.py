"""
Counter flow flame
 
Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""
 
import cantera as ct
import matplotlib.pyplot as plt
 
import numpy as np
 
 
# Input parameters
p = 101325  # pressure
tin_f = 293.0  # fuel inlet temperature
tin_o = 293.0  # oxidizer inlet temperature
mdot_o = 1.1996*0.1  # kg/m^2/s
mdot_f = 0.6655*0.1  # kg/m^2/s
 
comp_o = 'O2:0.21, N2:0.79'  # air composition
comp_f = 'CH4:1'  # fuel composition
 
 
width = 0.02  # Distance between inlets is 2 cm
 
loglevel = 1  # amount of diagnostic output (0 to 5)
 
# Create the gas object used to evaluate all thermodynamic, kinetic, and
# transport properties.
gas = ct.Solution('gri30.yaml')
gas.TP = gas.T, p
 
# Create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right.
f = ct.CounterflowDiffusionFlame(gas, width=width)
 
# Set the state of the two inlets
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f
 
f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o
 
# Set the boundary emissivities
f.boundary_emissivities = 0.0, 0.0
# Turn radiation off
f.radiation_enabled = False
 
#f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)
f.set_refine_criteria(ratio=2, slope=0.03, curve=0.08, prune=0.02)
 
# Solve the problem
f.solve(loglevel, auto=True)
 
 
# Calculate heat release rate
Q = -sum(f.standard_enthalpies_RT * f.net_production_rates, 0) * ct.gas_constant * f.T
# Save heat release rate
np.savetxt('Q.csv', list(zip(f.grid,Q)), delimiter=',')
 
# Save thermal conductivity
np.savetxt('kappa.csv', list(zip(f.grid,f.thermal_conductivity)), delimiter=',')
 
# Save viscosity
np.savetxt('mu.csv', list(zip(f.grid,f.viscosity)), delimiter=',')
 
# Save specific heat capacity
np.savetxt('cp.csv', list(zip(f.grid,f.cp_mass)), delimiter=',')
 
 
f.show_solution()
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('counteFlowFlame.h5', mode='w')
except ImportError:
    f.save('counteFlowFlame.yaml')
 
# write the velocity, temperature, and mass fractions to a CSV file
f.write_csv('counteFlowFlameResults.csv', species='Y', quiet=False)
 
f.show_stats(0)
 
# Plot Temperature
figTemperatureModifiedFlame = plt.figure()
plt.plot(f.flame.grid, f.T, label='Temperature')
plt.title('counteFlowFlame Temperature')
plt.ylim(0,2500)
plt.xlim(0.00, 0.02)
plt.legend()
plt.savefig('./counteFlowFlame.pdf')
