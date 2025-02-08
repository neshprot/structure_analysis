from reader import initialise
from octree import SpaceSeparation
from interactions import Config, Interactions
from utils import *
import os
import re


def get_file_names_from_folder(folder_path):
    """Возвращает список имен файлов в указанной папке."""
    try:
        file_names = os.listdir(folder_path)
        # Фильтруем, чтобы оставить только файлы (не каталоги)
        file_names = [f for f in file_names if os.path.isfile(os.path.join(folder_path, f))]
        return file_names
    except FileNotFoundError:
        print(f"Ошибка: Папка '{folder_path}' не найдена.")
        return []


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

def custom_sort_key(obj):
    """
    Разбивает строку объекта на последовательность букв и чисел для сортировки.
    """
    text = str(obj)
    parts = []
    for part in re.findall(r"(\d+|[^\d]+)", text):
        if part.isdigit():
            parts.append(int(part))
        else:
            parts.append(part)
    return parts

def complex_sort_key(key):
    parts = key[0].split(', ')
    first_part = parts[0].split(' - ')
    second_part = parts[1].split(' - ')
    return int(first_part[0]), int(second_part[0])

folder = "compare"
files = get_file_names_from_folder(folder)
pdb_names = sorted(files, key=custom_sort_key)

hbond_protein_dict = {}
electrostatic_protein_dict = {}
for pdb_name in pdb_names:
    protein = one_prot(f'compare/{pdb_name}')
    hbond_protein_dict, electrostatic_protein_dict = get_inteructions(protein, pdb_name, hbond_protein_dict, electrostatic_protein_dict)

sorted_items = sorted(hbond_protein_dict.items(), key=complex_sort_key)
create_excel(sorted_items, pdb_names, 'hbonds')

sorted_items = sorted(electrostatic_protein_dict.items(), key=complex_sort_key)
create_excel(sorted_items, pdb_names, 'electrostatic')
