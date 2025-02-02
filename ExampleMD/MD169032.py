from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


# load pdb structure and set force fields
pdb = PDBFile('LowEnergy169032.pdb')
forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')

# add solvent
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=1.0*nanometers, ionicStrength=0.15*molar)

# local energy minimization
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

# save results
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('output169032.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
with open('state_min169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
with open('system_min169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(simulation.system))


########################## NVT Equilibration ######################################
with open('system_min169032.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())
with open('state_min169032.xml', 'r') as f:
    state = XmlSerializer.deserialize(f.read())
pdb = PDBFile('output169032.pdb')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setState(state)
simulation.reporters.append(DCDReporter('NVT169032.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(50000)
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('NVT169032.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
with open('state_NVT169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
with open('system_NVT169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(simulation.system))


########################## NPT Equilibration ######################################
with open('system_NVT169032.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())
with open('state_NVT169032.xml', 'r') as f:
    state = XmlSerializer.deserialize(f.read())
pdb = PDBFile('NVT169032.pdb')
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setState(state)
simulation.reporters.append(DCDReporter('NPT169032.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(50000)
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('NPT169032.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
with open('state_NPT169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
with open('system_NPT169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(simulation.system))
    
    
########################## NPT MD Run ######################################
with open('system_NPT169032.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())
with open('state_NPT169032.xml', 'r') as f:
    state = XmlSerializer.deserialize(f.read())
pdb = PDBFile('NPT169032.pdb')
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setState(state)
simulation.reporters.append(DCDReporter('MD169032.dcd', 5000))
simulation.reporters.append(StateDataReporter(stdout, 5000, step=True, potentialEnergy=True, temperature=True))
simulation.step(500000000)
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('MD169032.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
with open('state_MD169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
with open('system_MD169032.xml', 'w') as f:
    f.write(XmlSerializer.serialize(simulation.system))

