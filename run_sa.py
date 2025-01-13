import time
from reader import initialise
import numpy as np
from octree import SpaceSeparation
from interactions import Config
import sys

start_time = time.time()
protein = initialise('ff.top', 'm1_equil.pdb')
atoms_in_prot = []
min_coord = np.array([float('inf'), float('inf'), float('inf')])
max_coord = np.array([-float('inf'), -float('inf'), -float('inf')])
# output of all class elements
for value in vars(protein).values():
    for k, v in value.items():
        for atom in v.atoms.values():
            if atom.name[0] in ['H', 'O', 'N']:
                atoms_in_prot.append(atom)
        coords = [atom.coordinates for atom in v.atoms.values()]
        min_coord = np.min(coords + [min_coord], axis=0)
        max_coord = np.max(coords + [max_coord], axis=0)

cubes = SpaceSeparation(atoms_in_prot, 5, min_coord, max_coord)

config = Config(1, 25, 30)

interaction = cubes.cube_cluster_search_action(protein.residues[219].atoms['HN'], 'hydrogen', config)

h_bonds = []

for res in protein.residues.values():
    for atom in res.atoms.values():
        interaction = cubes.cube_cluster_search_action(atom, 'hydrogen', config)

        h_bonds.append(sorted(interaction, key=lambda inter: [inter.atom1.name[0], -inter.squared_distance, inter.angle]))

print("time elapsed: {:.3f}s".format(time.time() - start_time))

