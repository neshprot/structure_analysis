import numpy as np
import math


class Config:
    def __init__(self, param1, param2, param3):
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
        if atom1.name[0] in ['H'] and atom2.name[0] in ['N', 'O', 'H']:
            point1 = atom1.coordinates
            donator = atom1.residue.bonds[atom1.name][0]
            point2 = atom1.residue.atoms[donator].coordinates
            point3 = atom2.coordinates
        elif atom2.name[0] in ['H'] and atom1.name[0] in ['N', 'O', 'H']:
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


class InteractionFactory:

    @staticmethod
    def create_interaction(interaction_type, *args, **kwargs):
        if interaction_type == "hydrogen":
            return HydrogenBond(*args, **kwargs)
