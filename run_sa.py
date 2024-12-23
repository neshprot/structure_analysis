import time
from reader import initialise
import numpy as np
from octree import SpaceSeparation

start_time = time.time()
protein = initialise('ff.top', 'm1_equil.pdb')
atoms_in_prot = []
min_coord = np.array([float('inf'), float('inf'), float('inf')])
max_coord = np.array([-float('inf'), -float('inf'), -float('inf')])
# output of all class elements
for value in vars(protein).values():
    for k, v in value.items():
        atoms_in_prot.extend(v.atoms.values())
        coords = [atom.coordinates for atom in v.atoms.values()]
        min_coord = np.min(coords + [min_coord], axis=0)
        max_coord = np.max(coords + [max_coord], axis=0)

cubes = SpaceSeparation(atoms_in_prot, 5, min_coord, max_coord)

s = 0
for atom in protein.residues[215].atoms.values():
    s += len(cubes.cube_cluster_search(atom, 1, 5))
print("time elapsed: {:.3f}s".format(time.time() - start_time))

