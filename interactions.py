from traceback import print_tb

import numpy as np
import math
from scipy.optimize import linear_sum_assignment


class Config:
    def __init__(self, param1, param2, param3):
        """

        :param param1: lower cutoff H-bond distance
        :param param2: upper cutoff H-bond distance
        :param param3: angle value cutoff
        """
        self.lower_cutoff = param1
        self.upper_cutoff = param2
        self.cutoff_angle = param3


class HydrogenBond:
    def __init__(self, atom1, atom2, config):
        self.atom1 = atom1
        self.atom2 = atom2

    def __new__(cls, atom1, atom2, config):
        lower_cutoff = config.lower_cutoff
        upper_cutoff = config.upper_cutoff
        if atom1.residue.res_id == atom2.residue.res_id:
            return None
        obj = super().__new__(cls)
        sq_dist = obj.calculate_squared_distance(atom1, atom2)
        if not lower_cutoff <= sq_dist <= upper_cutoff:
            return None
        angle = obj.calculate_angle(atom1, atom2)
        if not angle:
            return None
        if angle >= config.cutoff_angle:
            return None
        obj.squared_distance = sq_dist
        obj.angle = angle
        return obj

    @staticmethod
    def calculate_squared_distance(atom1, atom2):
        squared_distance = np.sum((atom1.coordinates - atom2.coordinates) ** 2)
        return squared_distance

    @staticmethod
    def calculate_angle(atom1, atom2):
        if len(atom1.residue.bonds[atom1.name]) > 1 or len(atom2.residue.bonds[atom2.name]) > 1:
            if atom1.residue.res_name != 'HOH' and atom2.residue.res_name != 'HOH':
                return None

        if atom1.name[0] in ['H'] and atom2.name[0] in ['N', 'O']:
            point1 = atom1.coordinates
            donator = atom1.residue.bonds[atom1.name][0]
            point2 = atom1.residue.atoms[donator].coordinates
            point3 = atom2.coordinates
        elif atom2.name[0] in ['H'] and atom1.name[0] in ['N', 'O']:
            point1 = atom2.coordinates
            donator = atom2.residue.bonds[atom2.name][0]
            point2 = atom2.residue.atoms[donator].coordinates
            point3 = atom1.coordinates
        else:
            return None

        vector1 = np.array(point1) - np.array(point2)
        vector2 = np.array(point3) - np.array(point2)

        # 2. Скалярное произведение
        dot_product = np.dot(vector1, vector2)

        # 3. Длины векторов
        mag_vector1 = np.linalg.norm(vector1)
        mag_vector2 = np.linalg.norm(vector2)

        # 4. Косинус угла
        cos_angle = dot_product / (mag_vector1 * mag_vector2)

        # 5. Угол (в радианах) и преобразование в градусы
        angle_rad = math.acos(cos_angle)
        angle_deg = math.degrees(angle_rad)
        return angle_deg


class ElectostaticBond:
    # atom1 - C in COO
    # atom2 - N in NH2
    def __init__(self, atom1, atom2, config):
        self.atom1 = atom1
        self.atom2 = atom2
    def __new__(cls, atom1, atom2, config):
        lower_cutoff = config.lower_cutoff
        upper_cutoff = config.upper_cutoff
        if atom1.residue.res_id == atom2.residue.res_id:
            return None
        if not atom2.name.startswith('N'):
            return None
        obj = super().__new__(cls)
        atom_o = list(filter(lambda s: s.startswith('O'), atom1.residue.bonds[atom1.name]))
        o_atoms = list(map(atom1.residue.atoms.get, atom_o, [None]*len(atom_o)))
        atom_h = list(filter(lambda s: s.startswith('H'), atom2.residue.bonds[atom2.name]))
        h_atoms = list(map(atom2.residue.atoms.get, atom_h, [None]*len(atom_h)))
        energy = 0
        for o in o_atoms:
            for h in h_atoms:
                sq_dist = obj.calculate_squared_distance(o, h)
                if not lower_cutoff <= sq_dist <= upper_cutoff:
                    return None
                energy += sq_dist
        obj.energy = energy
        return obj

    @staticmethod
    def calculate_squared_distance(atom1, atom2):
        squared_distance = np.sum((atom1.coordinates - atom2.coordinates) ** 2)
        return squared_distance


class InteractionFactory:

    @staticmethod
    def create_interaction(interaction_type, *args, **kwargs):
        if interaction_type == "hydrogen":
            return HydrogenBond(*args, **kwargs)
        if interaction_type == "electrostatic":
            return ElectostaticBond(*args, **kwargs)

class Interactions:
    @staticmethod
    def find_interaction_hydrogen(config, inter_name: str, separation, atoms, protein):
        def find_related_h_value(data, input_key):
            """
            Находит и выводит ключ-значение пару для ключа начинающегося на H,
            отличного от введенного ключа, если таковой имеется.
            :param data: Словарь, в котором ищем.
            :param input_key: Введенный ключ, начинающийся с H.
            :return: Вывод на консоль
            """
            h_keys = [key for key in data if key.startswith("H")]

            if not h_keys:
                return

            if input_key not in h_keys:
                return

            other_h_key = None
            for key in h_keys:
                if key != input_key:
                    other_h_key = key
                    break

            if other_h_key:
                return data[other_h_key]
            else:
                return

        inter_atoms = []
        for atom in atoms:
            nearest_atoms = separation.cube_cluster_search(atom, inter_name, config)
            for int_atom in nearest_atoms:
                if int_atom.name[0] in ['H']:
                    interaction_class = InteractionFactory.create_interaction(inter_name, atom, int_atom, config)
                    if interaction_class:
                        inter_atoms.append([interaction_class, atom.atom_id, int_atom.atom_id])

        edges = sorted(inter_atoms, key=lambda x: x[0].squared_distance)
        matched_nodes = set()
        result = []

        for weight, node1, node2 in edges:
            if node1 not in matched_nodes and node2 not in matched_nodes:
                result.append((weight, node1, node2))
                matched_nodes.add(node1)
                matched_nodes.add(node2)
                # only node1 cause node1 is donor and node2 acceptor
                if protein.atoms[node1].residue.res_name == 'HOH' and protein.atoms[node1].name[0] == 'H':
                    second_atom = find_related_h_value(protein.atoms[node1].residue.atoms,
                                                       protein.atoms[node1].name).atom_id
                    matched_nodes.add(second_atom)

        for inter, fisrt_a, second_a in result:
            if 'C' in inter.atom1.residue.bonds[inter.atom1.name] or 'N' in inter.atom1.residue.bonds[inter.atom1.name]:
                continue
            protein.atoms[fisrt_a].hbond = inter
            protein.atoms[second_a].hbond = inter

    @staticmethod
    def find_interaction_electrostatic(config, inter_name: str, separation, atoms, protein):
        def find_related_h_value(data, input_key):
            """
            Находит и выводит ключ-значение пару для ключа начинающегося на H,
            отличного от введенного ключа, если таковой имеется.
            :param data: Словарь, в котором ищем.
            :param input_key: Введенный ключ, начинающийся с H.
            :return: Вывод на консоль
            """
            h_keys = [key for key in data if key.startswith("H")]

            if not h_keys:
                return

            if input_key not in h_keys:
                return

            other_h_key = None
            for key in h_keys:
                if key != input_key:
                    other_h_key = key
                    break

            if other_h_key:
                return data[other_h_key]
            else:
                return

        inter_atoms = []
        for atom in atoms:
            nearest_atoms = separation.cube_cluster_search(atom, inter_name, config)
            for int_atom in nearest_atoms:
                if int_atom.name[0] in ['N']:
                    interaction_class = InteractionFactory.create_interaction(inter_name, atom, int_atom, config)
                    if interaction_class:
                        inter_atoms.append([interaction_class, atom.atom_id, int_atom.atom_id])

        edges = sorted(inter_atoms, key=lambda x: x[0].energy)
        matched_nodes = set()
        result = []

        for weight, node1, node2 in edges:
            if node1 not in matched_nodes and node2 not in matched_nodes:
                result.append((weight, node1, node2))
                matched_nodes.add(node1)
                matched_nodes.add(node2)

        for inter, fisrt_a, second_a in result:
            protein.atoms[fisrt_a].electrostatic = inter
            protein.atoms[second_a].electrostatic = inter
