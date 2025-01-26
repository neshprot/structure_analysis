import time
from reader import initialise
from octree import SpaceSeparation
from interactions import Config, Interactions
from utils import *


start_time = time.time()

def one_prot(file_name):
    protein = initialise('ff.top', file_name)

    atoms_in_hbonds, hbonds_donors, atoms_in_qq, qq_negative, min_coord, max_coord = atoms_for_hbond(protein)
    cubes_hbonds = SpaceSeparation(atoms_in_hbonds, 5, min_coord, max_coord)
    config_hbonds = Config(1, 25, 30)

    cubes_qq = SpaceSeparation(atoms_in_qq, 6, min_coord, max_coord)
    config_qq = Config(1, 36, 0)

    interactions = Interactions()
    interactions.find_interaction_hydrogen(config_hbonds, 'hydrogen', cubes_hbonds, hbonds_donors, protein)

    interactions.find_interaction_electrostatic(config_qq, 'electrostatic', cubes_qq, qq_negative, protein)
    return protein

pdb_names = []
i = 1
while i <= 400:
    pdb_names.append(f'5ZIM_b_{i}.pdb')
    i += 10

hbond_protein_dict = {}
for pdb_name in pdb_names:
    protein = one_prot(f'compare/{pdb_name}')
    hbond_protein_dict = get_inteructions(protein, pdb_name, hbond_protein_dict)


def complex_sort_key(key):
    parts = key[0].split(', ')
    first_part = parts[0].split(' - ')
    second_part = parts[1].split(' - ')
    return (int(first_part[0]), int(second_part[0]))

sorted_items = sorted(hbond_protein_dict.items(), key=complex_sort_key)
create_excel(sorted_items, pdb_names)

print("time elapsed: {:.3f}s".format(time.time() - start_time))
"""
for idx, res in sorted_items:
    res2 = protein2.residues[idx]
    for a, b in zip(res.atoms.values(), res2.atoms.values()):
        if a.hbond and (not b.hbond):
            print(f"Исчезла у residue {idx}: {format_interaction(a.hbond, idx)}")
            continue
        elif (not a.hbond) and b.hbond:
            print(f"Появилась у residue {idx}: {format_interaction(b.hbond, idx)}")
            continue
        elif a.hbond and b.hbond:
            if format_interaction(a.hbond, idx) != format_interaction(b.hbond, idx):
                print(f'Связь поменялась у residue {idx}: c {format_interaction(a.hbond, idx)} на {format_interaction(b.hbond, idx)}')

print("time elapsed: {:.3f}s".format(time.time() - start_time))

with open('logout', 'w') as ouf:
    sorted_items = sorted(protein.residues.items())
    for idx, res in sorted_items:
        if not any([atom.hbond for atom in res.atoms.values()]):
            continue
        ouf.write(f'Residue: {idx} - {res.res_name}\n')
        for a in res.atoms.values():
            if a.hbond:
                j = a.hbond
                ouf.write(f'{j.atom1.residue.res_id} - {j.atom1.name}, {j.atom2.residue.res_id} - {j.atom2.name}\n')
        ouf.write('\n')

with open('logout_qq', 'w') as ouf:
    sorted_items = sorted(protein.residues.items())
    for idx, res in sorted_items:
        if not any([atom.electrostatic for atom in res.atoms.values()]):
            continue
        ouf.write(f'Residue: {idx} - {res.res_name}\n')
        for a in res.atoms.values():
            if a.electrostatic:
                j = a.electrostatic
                ouf.write(f'{j.atom1.residue.res_id} - {j.atom1.name}, {j.atom2.residue.res_id} - {j.atom2.name}\n')
        ouf.write('\n')
"""
