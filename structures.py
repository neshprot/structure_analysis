from __future__ import annotations

import copy

import numpy as np
from numpy.typing import NDArray


class Atom:
    def __init__(self, residue: Residue):
        """
        Args:
            atom_name - name of atom
            coordinates - coordinates of atom
            :param residue - link to class residue, which consists this atom
        """
        self.__atom_id = None
        self.__name = None
        self.__coordinates = np.array([None, None, None])
        self.__residue = None
        self.__interaction = []

    @property
    def atom_id(self):
        return self.__atom_id

    @atom_id.setter
    def atom_id(self, value):
        self.__atom_id = value

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, atom_name: str) -> None:
        self.__name = atom_name

    @property
    def coordinates(self) -> NDArray[np.float64]:
        return self.__coordinates

    @coordinates.setter
    def coordinates(self, atom_coordinates: NDArray[np.float64]) -> None:
        self.__coordinates = atom_coordinates

    @property
    def residue(self) -> Residue:
        return self.__residue

    @residue.setter
    def residue(self, value):
        self.__residue = value

    @property
    def interaction(self) -> object:
        if self.__interaction:
            return self.__interaction[0]
        return None

    @interaction.setter
    def interaction(self, value: list):
        self.__interaction = value

    def back_interaction(self):
        return self.__interaction


class Residue:
    def __init__(self):
        """
        Args:
            res_name - name of residue
            bonds - graph with bonds (actually dictionary)
            res_id - residue id
            atoms - name of atoms associate with classes
        """
        self.__res_name = None
        self.__res_id = None
        self.__bonds = {}
        self.__atoms = {}

    @property
    def res_name(self) -> str:
        return self.__res_name

    @res_name.setter
    def res_name(self, name: str) -> None:
        self.__res_name = name

    @property
    def res_id(self) -> int:
        return self.__res_id

    @res_id.setter
    def res_id(self, num: int) -> None:
        self.__res_id = num

    @property
    def bonds(self) -> dict:
        return self.__bonds

    @bonds.setter
    def bonds(self, atom_pair: tuple[str, str]) -> None:
        self.__bonds.setdefault(atom_pair[0], []).append(atom_pair[1])
        self.__bonds.setdefault(atom_pair[1], []).append(atom_pair[0])

    @property
    def atoms(self) -> dict:
        return self.__atoms

    @atoms.setter
    def atoms(self, type_atom: tuple[str, Atom]) -> None:
        self.__atoms.update({type_atom[0]: type_atom[1]})

    def initialize_atoms(self, atom_name: str) -> None:
        """
        Initialize types of atoms in residue
        :param atom_name: name of atom
        :return: None
        """
        new_atom = Atom(self)
        new_atom.name = atom_name
        self.__atoms[atom_name] = new_atom
    '''
    def __copy__(self):
        new_res = Residue()
        new_res.__bonds = self.bonds
        new_res.__res_name = self.res_name
        new_res.__res_id = self.res_id
        for atom_type, atom in self.atoms.items():
            new_res.atoms = [atom_type, copy.deepcopy(atom)]
        return new_res
    '''

    def change_to_alt_res(self, alt_res: Residue) -> None:
        self.__res_name = alt_res.res_name
        for_adding = list(set(alt_res.atoms.keys()) - set(self.atoms.keys()))
        for_deleting = list(set(self.atoms.keys()) - set(alt_res.atoms.keys()))
        for atom_name in for_adding:
            self.initialize_atoms(atom_name)
        keys_to_delete = []
        bonds = copy.deepcopy(self.__bonds)
        for atom_name in for_deleting:
            del self.atoms[atom_name]
            for key, value in self.__bonds.items():
                if key == atom_name:
                    keys_to_delete.append(atom_name)
                else:
                    value[:] = [item for item in value if item != atom_name]
            for key in keys_to_delete:
                del bonds[key]
        for key, value in alt_res.__bonds.items():
            bonds.setdefault(key, []).extend(value)
        self.__bonds = bonds


class Protein:
    def __init__(self):
        """
        residues - all residues in protein
        """
        self.__residues = {}
        self.atoms = {}

    @property
    def residues(self):
        return self.__residues

    def add_residue(self, template_res: Residue, res_id: int) -> None:
        """
        add new residue to the protein, we should do it every time when check new residue id in pdb file
        :param template_res: name of residue
        :param res_id: it's id in protein
        :return: None
        """
        residue = copy.deepcopy(template_res)
        residue.res_id = res_id
        self.__residues[res_id] = residue
