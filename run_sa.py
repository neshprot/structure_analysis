import time
from reader import initialise
import numpy as np
from octree import SpaceSeparation
from interactions import Config, Interactions

start_time = time.time()
protein = initialise('ff.top', 'm1_equil.pdb')

def atoms_for_hbond(protein):
    def check_hbond(atom):
        if atom.name[0] in ['H', 'O', 'N']:
            if atom.residue.res_name != 'HOH':
                if len(atom.residue.bonds[atom.name]) != 1:
                    return False
            return True
        return False

    atoms_in_hbonds = []
    hbonds_donors = []
    atoms_in_qq = []
    qq_negative = []
    min_coord = np.array([float('inf'), float('inf'), float('inf')])
    max_coord = np.array([-float('inf'), -float('inf'), -float('inf')])
    for v in protein.residues.values():
        for atom in v.atoms.values():
            protein.atoms.update({atom.atom_id: atom})
            # skip backbone
            if 'C' in atom.residue.bonds[atom.name] or 'N' in atom.residue.bonds[atom.name]:
                continue
            # check hbonds
            if check_hbond(atom):
                atoms_in_hbonds.append(atom)
                if atom.name[0] in ['O', 'N']:
                    hbonds_donors.append(atom)

            # charge - charge
            if atom.residue.charge > 0:
                # check NH2-
                if atom.name[0] in ['N']:
                    hydrogens = list(filter(lambda s: s.startswith('H'), atom.residue.bonds[atom.name]))
                    if len(hydrogens) == 2:
                        atoms_in_qq.append(atom)
                        atoms_in_qq.extend(list(map(atom.residue.atoms.get, hydrogens, [None]*len(hydrogens))))
            elif atom.residue.charge < 0:
                if atom.name[0] in ['C']:
                    atom_o = list(filter(lambda s: s.startswith('O'), atom.residue.bonds[atom.name]))
                    if len(atom_o) == 2:
                        atoms_in_qq.extend(list(map(atom.residue.atoms.get, atom_o, [None]*len(atom_o))))
                        qq_negative.append(atom)
        coords = [atom.coordinates for atom in v.atoms.values()]
        min_coord = np.min(coords + [min_coord], axis=0)
        max_coord = np.max(coords + [max_coord], axis=0)
    return atoms_in_hbonds, hbonds_donors, atoms_in_qq, qq_negative, min_coord, max_coord

atoms_in_hbonds, hbonds_donors, atoms_in_qq, qq_negative, min_coord, max_coord = atoms_for_hbond(protein)
cubes_hbonds = SpaceSeparation(atoms_in_hbonds, 5, min_coord, max_coord)
config_hbonds = Config(1, 25, 30)

cubes_qq = SpaceSeparation(atoms_in_qq, 5, min_coord, max_coord)
config_qq = Config(1, 25, 0)

interactions = Interactions()
interactions.find_interaction_hydrogen(config_hbonds, 'hydrogen', cubes_hbonds, hbonds_donors, protein)

interactions.find_interaction_electrostatic(config_qq, 'electrostatic', cubes_qq, qq_negative, protein)

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
