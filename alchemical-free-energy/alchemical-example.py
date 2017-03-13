from simtk import openmm, unit

# Create a Lennard-Jones fluid
pressure = 80*unit.atmospheres
temperature = 120*unit.kelvin
collision_rate = 5/unit.picoseconds
timestep = 2.5*unit.femtoseconds
from openmmtools.testsystems import LennardJonesFluid
sigma = 3.4*unit.angstrom; epsilon = 0.238 * unit.kilocalories_per_mole
fluid = LennardJonesFluid(sigma=sigma, epsilon=epsilon)
[system, positions] = [fluid.system, fluid.positions]

# Add a barostat
barostat = openmm.MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# Retrieve the NonbondedForce
forces = { force.__class__.__name__ : force for force in system.getForces() }
nbforce = forces['NonbondedForce']

# Add a CustomNonbondedForce to handle only alchemically-modified interactions
alchemical_particles = set([0])
chemical_particles = set(range(system.getNumParticles())) - alchemical_particles
energy_function = 'lambda*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'
energy_function += 'reff_sterics = sigma*(0.5*(1.0-lambda) + (r/sigma)^6)^(1/6);'
energy_function += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
custom_force = openmm.CustomNonbondedForce(energy_function)
custom_force.addGlobalParameter('lambda', 1.0)
custom_force.addPerParticleParameter('sigma')
custom_force.addPerParticleParameter('epsilon')
for index in range(system.getNumParticles()):
    [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
    custom_force.addParticle([sigma, epsilon])
    if index in alchemical_particles:
        nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)
custom_force.addInteractionGroup(alchemical_particles, chemical_particles)
system.addForce(custom_force)

# Create a context
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
context = openmm.Context(system, integrator)
context.setPositions(positions)

# Minimize energy
print('Minimizing energy...')
openmm.LocalEnergyMinimizer.minimize(context)

# Collect data
nsteps = 2500 # number of steps per sample
niterations = 50 # number of samples to collect per alchemical state
import numpy as np
lambdas = np.linspace(1.0, 0.0, 10) # alchemical lambda schedule
nstates = len(lambdas)
u_kln = np.zeros([nstates,nstates,niterations], np.float64)
kT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * integrator.getTemperature()
for k in range(nstates):
    for iteration in range(niterations):
        print('state %5d iteration %5d / %5d' % (k, iteration, niterations))
        # Set alchemical state
        context.setParameter('lambda', lambdas[k])
        # Run some dynamics
        integrator.step(nsteps)
        # Compute energies at all alchemical states
        for l in range(nstates):
            context.setParameter('lambda', lambdas[l])
            u_kln[k,l,iteration] = context.getState(getEnergy=True).getPotentialEnergy() / kT

# Estimate free energy of Lennard-Jones particle insertion
from pymbar import MBAR, timeseries
# Subsample data to extract uncorrelated equilibrium timeseries
N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
for k in range(nstates):
    [nequil, g, Neff_max] = timeseries.detectEquilibration(u_kln[k,k,:])
    indices = timeseries.subsampleCorrelatedData(u_kln[k,k,:], g=g)
    N_k[k] = len(indices)
    u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T
# Compute free energy differences and statistical uncertainties
mbar = MBAR(u_kln, N_k)
[DeltaF_ij, dDeltaF_ij, Theta_ij] = mbar.getFreeEnergyDifferences()

print('DeltaF_ij (kT):')
print(DeltaF_ij)
print('dDeltaF_ij (kT):')
print(dDeltaF_ij)
