# Amber setup including Extra Particles - HSP90 with multisite-Mg2+ example

OpenMM allows users to model their systems using Amber and provide `prmtop`, and `inpcrd` files as input. This allows users more familiar with the Amber modeling environment to continue using their setup tools, while harnessing the speed and versatility of OpenMM. This also allows the use of non-standard force fields that have been published for use with Amber, such as metal ion models where dummy atoms are applied to mimic particular coordination geometries. Furthermore, this is facilitated by OpenMM's support for Extra Particles - i.e. particles that are not ordinary atoms, such as these metal dummy atoms, virtual sites in many water models etc.

This example illustrates the use of Amber's `tleap` to set up a simulation of the Heat Shock Protein 90 (HSP90, Uniprot: [`P07900`](http://www.uniprot.org/uniprot/P07900), including the ADP-Mg2+ complex, where we will model the Mg2+ cation using the Sept Lab octahedral multisite model ([`DOI:10.1021/ct400177g`](https://doi.org/10.1021/ct400177g), [`Sept Lab`](http://septlab.engin.umich.edu/multisite-ions.html)), and use the Carlson parameters for ADP ([`DOI:10.1002/jcc.10262`](https://doi.org/10.1002/jcc.10262), downloaded from [`Bryce group Amber parameter database`](http://research.bmh.manchester.ac.uk/bryce/amber)).

We begin from the 1BYQ PDB file, add missing residues (only those in the middle of the chain) and missing heavy atoms using PDBFixer. Using MDTraj, residue and atom naming of the Mg2+ ion and the ADP are fixed to match those in the parameter files, and dummy atoms of the multisite Mg2+ are added. `CONECT` bonds between Mg2+ and ADP are deleted. Crystallographic waters are preserved for coordination to the Mg2+. 

Finally, `tleap` is run to add hydrogens and parametrize. We use the `ff99SBildn` force field, `mag.lib` and `frcmod_mg.mmg` files from the multisite Mg2+ model (downloaded from [`Sept Lab`](http://septlab.engin.umich.edu/multisite-ions.html)), `ADP.prep` and `frcmod.phos` files for the ADP parameters (downloaded from [`Bryce group Amber parameter database`](http://research.bmh.manchester.ac.uk/bryce/amber)). The `prmtop` and `inpcrd` files are saved for simulation in OpenMM. The multisite Mg2+ model requires a correction to the Lennard-Jones B-coefficients in the `prmtop` file, this is done by a short Python script (obtained by email from Prof. David Sept, [`Sept Lab`](http://septlab.engin.umich.edu)).

```python
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
import mdtraj as md
import os
import numpy as np

# clean up the original PDB file and add missing residues and heavy atoms
fixer = PDBFixer('pdb1byq.ent')
fixer.findMissingResidues()

# only add missing residues in the middle of the chain, do not add terminal ones
chains = list(fixer.topology.chains())
keys = fixer.missingResidues.keys()
missingResidues = dict()
for key in keys:
    chain = chains[key[0]]
    if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
        missingResidues[key] = fixer.missingResidues[key]
fixer.missingResidues = missingResidues

fixer.findMissingAtoms()
fixer.addMissingAtoms()

PDBFile.writeFile(fixer.topology, fixer.positions, open('1byq_fixed.pdb', 'w'))

traj = md.load('1byq_fixed.pdb')

# fix residue names
# in adp - change ` to * in atom names
# add Mg dummy atoms
for residue in traj.top.residues:
    if residue.name == 'MG':
        residue.name = 'MMG'
        residue.atom(0).name = 'XZ'
        for i in range(6):
            traj.top.add_atom('D%s' % (i+1), md.element.Element.getBySymbol('Mg'), residue)
        # add dummy atoms
        central_index = residue.atom(0).index
        central_positions = traj.xyz[0,central_index]
        # positions from magnesium_ms.pdb, Sept Lab: http://septlab.engin.umich.edu/multisite-ions.html
        dummy_positions = central_positions + [
            [0.0, 0.0, 0.09],
            [0.09, 0.0, 0.0],
            [0.0, 0.0, -0.09],
            [-0.09, 0.0, 0.0],
            [0.0, -0.09, 0.0],
            [0.0, 0.09, 0.0]
        ]    
        traj.xyz = np.array([np.insert(traj.xyz[0], central_index+1, dummy_positions, axis=0)])    
    elif residue.name == 'ADP':
        residue.name = 'adp'
        for atom in residue.atoms:
            if atom.name[-1] == "'":
                atom.name = atom.name[:-1] + "*"
                
# remove adp - mg CONECT bonds
bonds = []

for bond in traj.top._bonds:
    if not (bond[0].residue.name == 'MMG' or bond[1].residue.name == 'MMG'):
        bonds.append(bond)
        
traj.top._bonds = bonds

traj.save('1byq_fixed2.pdb')

# save the tleap script to file
with open('leaprc.hsp90', 'w') as f:
    f.write('''
source oldff/leaprc.ff99SBildn
loadOff mag.lib
loadamberparams frcmod_mg.mmg
loadAmberPrep ADP.prep
loadamberparams frcmod.phos
x = loadPdb 1byq_fixed2.pdb
addIons x Na+ 0
solvateBox x TIP3PBOX 10.0
savePdb x topology.pdb
saveAmberParm x input.prmtop input.inpcrd
quit
''')

# run tleap
os.system('tleap -f leaprc.hsp90')
```

```
> python zeroBvalues.py 

This is a simple script designed to zero out the Lennard-Jones
B values for dummy atoms in our multisite ion model.

Please input topology file: input.prmtop
Modified parameter file will be input.prmtop.mod
Topology file has 49700 atoms and 22 atom types
D1 atom at position 3358 in atom list
D1 atom has type 16
Making substitution in LENNARD_JONES_BCOEF at positions: [121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 152, 169, 187, 206, 226, 247]
```

We load the `prmtop` and `inpcrd` files in, by creating `AmberPrmtopFile` and `AmberInpcrdFile` objects. Next, the `System` is created by calling the `createSystem` method on the `AmberPrmtopFile` object. Next, the `LangevinIntegrator` and the `Simulation` are set up, using the topology from the `AmberPrmtopFile` and positions from the `AmberInpcrdFile`. In this example we will use the `CUDA` platform, with mixed precision. The simulation is energy minimized and equilibrated for 100 steps. Reporters are attached and the production simulation propagated for 50 ns. 

```python
from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

# load in Amber input files
prmtop = app.AmberPrmtopFile('input.prmtop.mod')
inpcrd = app.AmberInpcrdFile('input.inpcrd')

# prepare system and integrator
system = prmtop.createSystem(nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

# prepare simulation
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(prmtop.topology, system, integrator, platform, 
    properties)
simulation.context.setPositions(inpcrd.positions)

# minimize
print('Minimizing...')
simulation.minimizeEnergy()

# equilibrate for 100 steps
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)

# append reporters
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=25000000, separator='\t'))

# run 50 ns of production simulation
print('Running Production...')
simulation.step(25000000)
print('Done!')
```
