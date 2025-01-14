from structures import Atom
from interactions import InteractionFactory, Config


class SpaceSeparation:
    def __init__(self, atoms: list[Atom], cell_size: int, min_coord, max_coord):
        """
        Args:
            cubes - list of sectors after partition
        :param atoms: all atoms in protein
        :param cell_size: size of one cell
        :param min_coord: min coordinates of atoms
        :param max_coord: max coordinates of atoms
        """
        self.__cell_size = cell_size
        self.__atoms = atoms
        self.__min_coord = min_coord
        self.__max_coord = max_coord
        self.__cubes = {}
        self.cube_comparison()

    def cube_comparison(self):
        for atom in self.__atoms:
            relative_position = [coord - start for coord, start in zip(atom.coordinates, self.__min_coord)]
            sector_index = [str(int(rel_p // 5)) for rel_p in relative_position]
            sec_id = "-".join(sector_index)
            self.__cubes.setdefault(sec_id, []).append(atom)

    def cube_cluster_search(self, core_atom: Atom, interaction_type: str, config: Config) -> list:
        sectors_around = [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1),
                          (-1, 0, 0), (-1, 0, 1), (-1, 1, -1), (-1, 1, 0),
                          (-1, 1, 1), (0, -1, -1), (0, -1, 0), (0, -1, 1),
                          (0, 0, -1), (0, 0, 1), (0, 1, -1),
                          (0, 1, 0), (0, 1, 1), (1, -1, -1), (1, -1, 0),
                          (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1),
                          (1, 1, -1), (1, 1, 0), (1, 1, 1), (0, 0, 0)]

        relative_position = [coord - start for coord, start in zip(core_atom.coordinates, self.__min_coord)]
        core_sector_index = [str(int(rel_p // self.__cell_size)) for rel_p in relative_position]
        """
        # sectors around central one
        for shift in sectors_around:
            sector_index = [str(int(a) - b) for a, b in zip(core_sector_index, shift)]
            sec_id = "-".join(sector_index)
            for atom in self.__cubes.get(sec_id, []):
                if atom.interaction:
                    continue
                interaction_class = InteractionFactory.create_interaction(interaction_type, core_atom, atom, config)
                if interaction_class:
                    interactions.append(interaction_class)
                    core_atom.interaction = interaction_class
                    atom.interaction = interaction_class
        return interactions
        """
        nearest_atoms = []
        for shift in sectors_around:
            sector_index = [str(int(a) - b) for a, b in zip(core_sector_index, shift)]
            sec_id = "-".join(sector_index)
            for atom in self.__cubes.get(sec_id, []):
                nearest_atoms.append(atom)
        return nearest_atoms
