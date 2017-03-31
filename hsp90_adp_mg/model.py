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
