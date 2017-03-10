# Alchemical free energy example

OpenMM's custom forces---which allow the programmer to express a potential algebraically, potentially with multiple parameters that can be adjusted on the fly---allow a great deal of flexibility and simplicity in encoding potentials while still achieving high performance on GPUs.
One common use of this facility is to convert standard interactions (such as Lennard-Jones potentials) into alchemically-modified potentials for the purposes of computing free energy differences.
The alchemical free energy code [YANK](http://github.com/choderalab/yank), for example, uses a variety of custom forces to represent alchemically-modified potentials for [protein-ligand alchemical binding free energy calculations](http://dx.doi.org/10.1007/s10822-013-9689-8).

As a simple example of how this is facilitated by custom forces, consider computing the chemical potential of liquid argon by estimating the free energy of inserting a Lennard-Jones particle 

```python
from si
# Create a Lennard-Jones fluid
system = openmm

```
