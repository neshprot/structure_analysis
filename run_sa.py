import time
from reader import initialise
import numpy as np
from octree import SpaceSeparation
from interactions import Config, InteractionFactory

start_time = time.time()
protein = initialise('ff.top', 'm1_equil.pdb')

atoms_in_prot = []
min_coord = np.array([float('inf'), float('inf'), float('inf')])
max_coord = np.array([-float('inf'), -float('inf'), -float('inf')])
for k, v in protein.residues.items():
    for atom in v.atoms.values():
        if atom.name[0] in ['H', 'O', 'N'] and len(atom.residue.bonds[atom.name]) == 1:
            atoms_in_prot.append(atom)
    coords = [atom.coordinates for atom in v.atoms.values()]
    min_coord = np.min(coords + [min_coord], axis=0)
    max_coord = np.max(coords + [max_coord], axis=0)
cubes = SpaceSeparation(atoms_in_prot, 5, min_coord, max_coord)
config = Config(1, 25, 30)


def find_interactions(config, inter_name: str, separation: SpaceSeparation, atoms):
    hbonds = {}
    for atom in atoms:
        interactions = []
        nearest_atoms = cubes.cube_cluster_search(atom, inter_name, config)
        for int_atom in nearest_atoms:
            if atom.interaction:
                continue
            interaction_class = InteractionFactory.create_interaction(inter_name, atom, int_atom, config)
            if interaction_class:
                interactions.append(interaction_class)
                int_atom.interaction = interaction_class
                atom.interaction = interaction_class
        if interactions:
            intr = sorted(interactions, key=lambda inter: [inter.atom1.name[0], inter.squared_distance, inter.angle])[0]
            hbonds.setdefault(intr.atom1.residue.res_id, {}).update({intr.atom1.name: intr})
            hbonds.setdefault(intr.atom2.residue.res_id, {}).update({intr.atom2.name: intr})
    return hbonds

hbonds = find_interactions(config, 'hydrogen', cubes, atoms_in_prot)

i = hbonds[215]
for n, j in i.items():
    print(f'{j.atom1.residue.res_id} - {j.atom1.name}, {j.atom2.residue.res_id} - {j.atom2.name}')
print("time elapsed: {:.3f}s".format(time.time() - start_time))
