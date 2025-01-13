import numpy as np
from structures import Residue, Protein
import copy


class Topology:
    def __init__(self):
        self.__residues = {}
        self.__association = {}

    @property
    def residues(self) -> dict:
        return self.__residues

    @residues.setter
    def residues(self, name_res: tuple[str, Residue]) -> None:
        self.__residues.update({name_res[0]: name_res[1]})

    @property
    def association(self):
        return self.__association

    def read_topology(self, top_file: str) -> None:
        with open(top_file, 'r') as inf:
            res_type = None
            for line in inf:
                line = line.strip()     # delete spaces in the beginning and the end
                # create new residues
                if line.startswith("RESI") or line.startswith('PRES'):
                    res_name = line.split()[1]
                    residue = Residue()
                    residue.res_name = res_name
                    self.residues = [res_name, residue]
                    res_type = res_name
                    continue
                # add atoms
                if line.startswith("ATOM"):
                    atom_line = line.split()
                    atom_name = atom_line[1]
                    self.residues[res_type].initialize_atoms(atom_name)
                    continue
                # add bonds
                if line.startswith("BOND"):
                    bonds_line = line.split('!')[0].split()[1::]
                    for atom1, atom2 in zip(bonds_line[::2], bonds_line[1::2]):
                        self.residues[res_type].bonds = [atom1, atom2]
                    continue
                if line.startswith("PATCH"):
                    patch_line = line.split()[2::2]
                    self.__association.setdefault(res_type, []).extend(patch_line)


class PdbReader:
    def __init__(self, topology, file_name):
        self.__topology = topology
        self.__protein = Protein()
        self.read_pdb(file_name)

    @property
    def topology(self):
        return self.__topology

    @property
    def protein(self):
        return self.__protein

    def read_pdb(self, file_name) -> None:
        id_type = None  # residue id
        with open(file_name, 'r') as inf:
            for line in inf:
                line = line.strip()  # delete spaces in the beginning and the end
                if line.startswith('ATOM'):
                    atom_line = line.split()[1:-3:]
                    res_id = int(atom_line[4])
                    atom_name = atom_line[1]
                    x_coord = float(atom_line[5])
                    y_coord = float(atom_line[6])
                    z_coord = float(atom_line[7])
                    cartesian_coords = [x_coord, y_coord, z_coord]
                    if res_id != id_type:
                        id_type = res_id
                        name_type = atom_line[2]
                        self.__protein.add_residue(self.__topology.residues[name_type], res_id)
                    if atom_name not in self.__protein.residues[id_type].atoms:
                        for alt_name in self.__topology.association[name_type]:
                            if atom_name in self.__topology.residues[alt_name].atoms:
                                name_type = alt_name
                                self.__protein.residues[id_type].change_to_alt_res(self.__topology.residues[name_type])
                                break
                        self.__protein.residues[id_type].res_id = res_id
                    self.__protein.residues[id_type].atoms[atom_name].coordinates = np.array(cartesian_coords)
                    self.__protein.residues[id_type].atoms[atom_name].residue = self.__protein.residues[id_type]


def initialise(top_file, pdb_file):
    protein_topology = Topology()
    protein_topology.read_topology(top_file)
    protein = PdbReader(protein_topology, pdb_file).protein
    return protein
