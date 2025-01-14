import time
from reader import initialise
import numpy as np
from octree import SpaceSeparation
from interactions import Config, Interactions

start_time = time.time()
protein = initialise('ff.top', 'm1_equil.pdb')

def atoms_for_hbond(protein):
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
    return atoms_in_prot, min_coord, max_coord

atoms_in_prot, min_coord, max_coord = atoms_for_hbond(protein)
cubes = SpaceSeparation(atoms_in_prot, 5, min_coord, max_coord)
config = Config(1, 25, 30)
interactions = Interactions()
interactions.find_interactions(config, 'hydrogen', cubes, atoms_in_prot)

print("time elapsed: {:.3f}s".format(time.time() - start_time))

with open('logout', 'w') as ouf:
    sorted_items = sorted(protein.residues.items())
    for idx, res in sorted_items:
        if not any([atom.interaction for atom in res.atoms.values()]):
            continue
        ouf.write(f'Residue: {idx} - {res.res_name}\n')
        for a in res.atoms.values():
            if a.interaction:
                j = a.interaction
                ouf.write(f'{j.atom1.residue.res_id} - {j.atom1.name}, {j.atom2.residue.res_id} - {j.atom2.name}\n')
        ouf.write('\n')
