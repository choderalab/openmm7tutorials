Using CHARMM input files in OpenMM and a CHARMM-GUI example
-----------------------------------------------------------

OpenMM can directly read CHARMM input files through the use of the `simtk.openmm.app` layer. This enables the use of all 
the powerful setup tools in the CHARMM ecosystem that a user might be familiar with such as the [CHARMM-GUI](http://onlinelibrary.wiley.com/doi/10.1002/jcc.20945/abstract), [VMD](http://www.sciencedirect.com/science/article/pii/0263785596000185?via%3Dihub), [CGenFF program](https://cgenff.paramchem.org/)
(DOI:[10.1002/jcc.21367](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21367/abstract))etc. This allows users who are 
already working in the CHARMM environment to harness the GPU speeds that OpenMM provides without having to modify their 
simulation system description files. 


OpenMM can also read CHARMM force field files. Therefore it is possible to use force fields that aren't already included 
in OpenMM such as the general CHARMM force field (CGenFF). For example, a user can generate an `str` file with the CGenFF program
for a ligand and load it into OpenMM. However, when using this feature, the `CharmmParameterSet` class needs to be used to 
load all the other CHARMM force field files as demonstrated in the example in listing 1. 


The example demonstrates how to use CHARMMM files that were generated with the CHARMM-GUI in an OpenMM script.The OpenMM app layer includes several
classes to load CHARMM files. The `CharmmPsfFile` class reads the `psf` file and instantiates a chemical structure on 
which one can then call the `createSystem()` method to creates an OpenMM system object. For the atomic coordinates, a 
regular `pdb` file can be used or the `CharmmCrdFile` or `CharmmRstFile` classes can be used to read CHARMM coordinate 
files. Files containing force field definitions come in a variety of formats such as `prm`, `par`, `top`, `rtf`, `inp` 
and `str`. These files are loaded into a `CharmmParameterSet` object which is then included as the first parameter when 
`createSystem()` is called on the chemical structure. For this example, the [membrane builder](http://dx.doi.org/10.1371/journal.pone.0000880) in the CHARMM-GUI was used
to generate the input files for the B2AR in a POPC lipid membrane. The membrane builder provides a straightforward way 
to go from the RCSB X-ray structure to the protein embedded in a membrane with all the relevant CHARMM input 
files.

```python
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

# Load CHARMM files
psf = CharmmPsfFile('step5_charmm2omm.psf')
pdb = PDBFile('step5_charmm2omm.pdb')
params = CharmmPsfFile('par_all36_prot.rtf', 'top_all36_prot.prm',
                       'par_all36_lipid.rtf', 'top_all36_lipid.prm',
                       'toppar_water_ion.str')

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

The CHARMM-GUI also generates a more elaborate set of OpenMM scripts to run equilibration and the production simulation
that are very straightforward to use. A tutorial that walks through this process is provided [here](https://github.com/ChayaSt/openmm7tutorials/tree/master/b2ar_membrane). 
When OpenMM is selected in the last step, CHARMM-GUI provides all the 
relevant openmm scripts in the downloaded tarball. The parameters and 
arguments for openmm objects and functions are provided in the `inp` 
files for all equilibration steps and production. This makes it simple 
to change parameters such as the time-step or electrostatic cut-off method. 