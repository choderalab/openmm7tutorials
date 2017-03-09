Using OpenMM7 to run a simulation of b2ar in a lipid bilayer.
=============================================================


CHARMM-GUI to generate the lipid bilayer.
-----------------------------------------

The pdb file used in this tutorial, `apo_snapshot.pdb` was taken from the reference structures for [Kohloff et.al. 2014](10.1038/NCHEM.1821), 
Provenance [3POG](http://www.rcsb.org/pdb/explore.do?structureId=3p0g).


1. Go to [CHARMM-GUI's membrane builder](http://www.charmm-gui.org/?doc=input/membrane) and upload the pdb file 
`apo_snapshot.pdb`. Select PDB format and click `Next Step`.
2. Stick with the default selection (Protein) and click `Next Step`. 
3. Select `Terminal group patching` and select termin `NTER` and `CTER`. 
4. Select `Disulfide bonds` for the pairs `184-190` and `106-191` and click `Next Step`. 
5. Select `Align the First Principal Axis Along Z` (the suggestion for a homo-oligomer) and click `Next Step`. 
6. Keep the default selections for `Heterogeneous Lipid`, `Rectangular` Box Type, `17.5 water thickness`. Select 
`Number of lipid components` and select `POPC` lipids with `55` lipids on the upperleaflet and lowerleaflet. Clike
`Next Step`. 
7. Keep all default selections (`Replacement method`, `check lipid ring (and protein surface) penetration`, 
`Include Ions 0.15 KCl` and `Ion Placing Method: Monte-Carlo`. Click `Next Step`. 
8. Click `Next Step` to assemble components. When you arrive at `Step 5`, select `OpenMM` as the Input generation
option. Keep all the other default options (`Generate grid information for PME FFT automatically` and `NPT ensemble`. 
Click `Next Step`. 

Equilibrating and running the simulation with OpenMM7.
------------------------------------------------------

From the `charmm-gui.tgz` extract `toppar/` directory and the `openmm/` directory. The `README` file in the `openmm/` 
directory is a csh script to run equilibration and the production simulation.

```
sch README
```


