from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

from openmmtools import states, mcmc, multistate
from openmmtools.states import ThermodynamicState, SamplerState

name = sys.argv[1]

# Load the NVT equilibrated system
with open('system_NVT'+name+'.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())
with open('state_NVT'+name+'.xml', 'r') as f:
    state = XmlSerializer.deserialize(f.read())

temperatures = [300.00, 304.48, 309.01, 313.60, 318.24, 322.94, 327.70, 332.52, 337.40, 342.34, 347.34, 352.43, 357.56, 362.75, 368.00, 373.32, 378.70, 384.15, 389.67, 395.26, 400.92, 406.65, 412.45, 418.32, 424.26, 430.28, 436.39, 442.56, 448.81, 455.13, 461.54, 468.01, 474.58, 481.22, 487.95, 494.79, 500.00] * kelvin
n_replicas = len(temperatures)
n_iterations = 50000

# Set the states, NPT ensemble
thermodynamic_states = [ThermodynamicState(system=system, temperature=T, pressure=1*bar) for T in temperatures]
sampler_state = SamplerState(positions=state.getPositions(), velocities=state.getVelocities(), box_vectors=state.getPeriodicBoxVectors())

# Use Langevin dynamics integrator
move = mcmc.LangevinSplittingDynamicsMove(timestep=2*femtosecond, collision_rate=1/picosecond, n_steps=1000)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, replica_mixing_scheme='swap-neighbors', number_of_iterations=n_iterations)

# Run the simulation and output the results in .nc file
reporter = multistate.MultiStateReporter('REMD'+name+'.nc', checkpoint_interval=1)
simulation.create(thermodynamic_states,sampler_state,storage=reporter)
simulation.run()
