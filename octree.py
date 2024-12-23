from structures import Atom
import numpy as np


class SpaceSeparation:
    def __init__(self, atoms: list[Atom], cell_size: int, min_coord, max_coord):
        self.__cell_size = cell_size
        self.__atoms = atoms
        self.__min_coord = min_coord
        self.__max_coord = max_coord
        self.__cubes = {}
        self.cube_comparison()

    def pairwise_comparison(self, core_atom, lower_cutoff, upper_cutoff):
        squared_lower_cutoff = lower_cutoff ** 2
        squared_upper_cutoff = upper_cutoff ** 2
        core_coord = core_atom.coordinates
        interaction_atoms = []
        for atom in self.__atoms:
            squared_distance = np.sum((core_coord - atom.coordinates) ** 2)
            if squared_lower_cutoff <= squared_distance <= squared_upper_cutoff:
                interaction_atoms.append(atom)
        return interaction_atoms

    def cube_comparison(self):
        for atom in self.__atoms:
            relative_position = [coord - start for coord, start in zip(atom.coordinates, self.__min_coord)]
            sector_index = [str(int(rel_p // 5)) for rel_p in relative_position]
            sec_id = "-".join(sector_index)
            self.__cubes.setdefault(sec_id, []).append(atom)

    def cube_cluster_search(self, core_atom, lower_cutoff, upper_cutoff):
        squared_lower_cutoff = lower_cutoff ** 2
        squared_upper_cutoff = upper_cutoff**2

        def check_interaction(core_c, p_atom):
            squared_distance = np.sum((core_c - p_atom.coordinates) ** 2)
            if squared_lower_cutoff <= squared_distance <= squared_upper_cutoff:
                return True
            return False

        sectors_around = [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1),
                          (-1, 0, 0), (-1, 0, 1), (-1, 1, -1), (-1, 1, 0),
                          (-1, 1, 1), (0, -1, -1), (0, -1, 0), (0, -1, 1),
                          (0, 0, -1), (0, 0, 1), (0, 1, -1),
                          (0, 1, 0), (0, 1, 1), (1, -1, -1), (1, -1, 0),
                          (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1),
                          (1, 1, -1), (1, 1, 0), (1, 1, 1)]

        core_coord = core_atom.coordinates
        interaction_atoms = []

        relative_position = [coord - start for coord, start in zip(core_atom.coordinates, self.__min_coord)]
        core_sector_index = [str(int(rel_p // self.__cell_size)) for rel_p in relative_position]
        core_sec_id = "-".join(core_sector_index)
        # central sector
        for atom in self.__cubes[core_sec_id]:
            if check_interaction(core_coord, atom):
                interaction_atoms.append(atom)
        # sectors around central one
        for shift in sectors_around:
            sector_index = [str(int(a) - b) for a, b in zip(core_sector_index, shift)]
            sec_id = "-".join(sector_index)
            for atom in self.__cubes.get(sec_id, []):
                if check_interaction(core_coord, atom):
                    interaction_atoms.append(atom)

        return interaction_atoms
