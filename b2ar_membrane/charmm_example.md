OpenMM can directly read CHARMM input files through the use of the `simtk.openmm.app` layer.
This enables the use of all the powerful setup tools in the CHARMM ecosystem such as the CHARMM-GUI,
VMD, CGenFF, etc. Since OpenMM can also read CHARMM force field files, it is possible to use force fields that
aren't already included in OpenMM.
The example in listing 1 demonstrates how to use CHARMMM files in an OpenMM sctipt.
The OpenMM app layer includes several classes to load CHARMM files. The `CharmmPsfFile` class reads the `psf` file and 
instantiates a chemical structure which is then used to create an OpenMM system. For the atomic coordinates, a regular
`pdb` file can be used or the `CharmmCrdFile` or `CharmmRstFile` classes can be used to read CHARMM coordinate files. 
Files containing force field definitions come in a variety of formats such as `prm`, `par`, `top`, `rtf`, `inp` and
`str`. These files are loaded into a `CharmmParameterSet` object which is then included as the first parameter when 
`createSystem()` is called on the chemical structure.

```
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

# Load CHARMM files
psf = CharmmPsfFile('input.psf')
pdb = PDBFile('input.pdb')
params = CharmmPsfFile('charmm36.rtf', 'charmm36.prm')

# Create an openmm system by calling createSystem on psf
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
         nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin,   # Temperature of head bath
                                1/picosecond, # Friction coefficient
                                0.002*picoseconds) # Time step
                     
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

# Set up the reporters to report energies every 1000 steps.
simulaiton.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                            potentialEnergy=True, temperature=True))
# run simulation
simulation.step(10000)

```

The CHARMM-GUI can generate a more elaborate set of OpenMM scripts to run 
simulations and they are very straightforward to use. As an example, consider
setting up a simulation of the B2AR GPCR in a membrane for simulation. 
To do this, we will us the membrane builder in the CHARMM-GUI. 
The membrane builder provides a straightforward way to go from the RCSB X-ray 
structure to the protein embedded in a membrane with all the relevant CHARMM input 
files. A tutorial that walks through this process is provided [here](https://github.com/ChayaSt/openmm7tutorials/tree/master/b2ar_membrane). 
When OpenMM is selected in the last step, CHARMM-GUI provides all the 
relevant openmm scripts in the downloaded tarball. The parameters and 
arguments for openmm objects and functions are provided in the `inp` 
files for all equilibration steps and production. This makes it simple 
to change parameters such as the time-step or electrostatic cut-off method. 