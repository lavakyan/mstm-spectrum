import yaml
import numpy as np
from mstm_studio.mstm_spectrum import Material
import zipfile
import mstm_studio_support as sup
import os
import matplotlib
matplotlib.use('TkAgg')
from tkinter import filedialog, messagebox, Menu
import json
import hashlib

global extract_dir
extract_dir = None      # Global variable for the extraction directory
valid_files = []        # Global list for storing validated files

def import_from_file_plot(filename):
    # Load YAML data from the file and process it
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)
    data_points = data['DATA'][0]['data'].strip().splitlines()
    wls = []
    n_values = []
    k_values = []

    # Process the data points
    for line in data_points:
        wl, n, k = map(float, line.split())
        wls.append(wl * 1000)  # Convert wavelength to nm
        n_values.append(n)
        k_values.append(k)
    wls = np.array(wls)
    n_values = np.array(n_values)
    k_values = np.array(k_values)
    
    parts = filename.split('\\')
    element = parts[-2]
    compound = parts[-1].split('.')[0] 
    name = f"{element}({compound})"

    return Material(
        file_name=name,
        wls=wls,
        nk=n_values + 1j * k_values 
    )

def get_file_hash(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def btLoadDatabaseClick(master=None):
    global extract_dir, valid_files
    valid_files = []
    extract_dir = None
    zip_file_path = filedialog.askopenfilename(title="Select a ZIP file", filetypes=[("ZIP files", "*.zip")])
    
    if not zip_file_path or not zip_file_path.lower().endswith('.zip'):
        messagebox.showerror("Error", "Please select a ZIP file.")
        return
    db_name = os.path.splitext(os.path.basename(zip_file_path))[0]
    cache_dir = "database_cache"
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    cache_file = os.path.join(cache_dir, f"{db_name}.json")
    if os.path.exists(cache_file):
        try:
            with open(cache_file, "r") as cache:
                valid_files = json.load(cache)
            messagebox.showinfo("Cache Loaded", "Cache loaded")
            extract_dir = os.path.dirname(cache_file)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load cache: {e}")
            return
    if not valid_files:
        extract_dir = os.path.join(os.getcwd(), db_name)
        try:
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to extract the file: {e}")
            return

        print(f"Material database extracted to {extract_dir}")
        def filter_valid_files(current_path):
            try:
                for root, _, files in os.walk(current_path):
                    for file in files:
                        if file.lower().endswith('.yml'):
                            item_path = os.path.join(root, file)
                            try:
                                mat = import_from_file_plot(item_path)
                                mat.get_n(250)
                                mat.get_k(250)
                                mat.get_n(1000)
                                mat.get_k(1000)
                                valid_files.append(item_path)
                            except Exception as e:
                                continue
            except Exception as e:
                messagebox.showerror("Error", f"Error while filtering files in {current_path}: {e}")

        filter_valid_files(extract_dir)

    try:
        with open(cache_file, "w") as cache:
            json.dump(valid_files, cache)
        print(f"Cache saved to {cache_file}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to save cache: {e}")

    messagebox.showinfo("Success!", "Database loaded successfully.")

def btLoadData(root, menu_button):
    global valid_files
    if extract_dir is None or not valid_files:
        messagebox.showerror('Warning', 'Database not imported')
    else:
        menu = Menu(root, tearoff=0)

        def build_menu(menu, valid_files):
            file_structure = {}
            for file_path in valid_files:
                relative_path = os.path.relpath(file_path, extract_dir)
                parts = relative_path.split(os.sep)
                current_dict = file_structure
                for part in parts[:-1]:
                    current_dict = current_dict.setdefault(part, {})
                current_dict[parts[-1]] = file_path

            def add_menu_items(menu, structure):
                for name, item in structure.items():
                    if isinstance(item, dict):
                        submenu = Menu(menu, tearoff=0)
                        if add_menu_items(submenu, item):
                            menu.add_cascade(label=name, menu=submenu)
                    else:
                        menu.add_command(label=name, command=lambda p=item: select_yml_file(p))
                return bool(menu.index("end") is not None)

            add_menu_items(menu, file_structure)

        def select_yml_file(file_path):
            mat = import_from_file_plot(file_path)
            key = sup.gen_mat_key()
            sup.add_material(key, mat)
            sup.update_materials_tree()

        build_menu(menu, valid_files)
        x = menu_button.winfo_rootx() + menu_button.winfo_width()
        y = menu_button.winfo_rooty()
        menu.post(x, y)
