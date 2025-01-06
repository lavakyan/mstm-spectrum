
# ~ import os
import sys
import numpy as np
from scipy import interpolate
from mstm_studio.mstm_spectrum import Material

try:
    from pathlib import Path
    # ~ import json
    # ~ import hashlib
    import zipfile
    import yaml
except:
    print('WARNING: reading from refractiveindex.info is disabled')
    print('required python modules: zipfile, yaml')

# ~ from tkinter import filedialog, messagebox, Menu
# ~ import mstm_studio.mstm_studio_support as sup


class RiiMaterial(Material):

    def __init__(self, archive_filename=''):
        '''
        Setup material from RefractiveIndexInfo database dump
        available online at url:
        <https://refractiveindex.info/download/database/>

        archive_filename: string
            path to the downloaded zip file of db dump

        Example of usage:
        >>> riimat = RiiMaterial('rii-database-2024-08-14.zip')
        >>> riimat.select('main', 'Ag', 'Babar')
        '''
        # ~ print(archive_filename)
        self.rii_db_items = dict()
        if archive_filename == '':
            # use system-dependedt default path
            if sys.platform == 'win32':
                appdata_path = Path.home() / 'AppData' / 'Roaming' / 'mstm_studio'
            else:
                appdata_path = Path.home() / '.mstm_studio'
            # find in appdata_path
            filenames = list(appdata_path.glob('rii-database-*.zip'))
            if len(filenames) > 0:
                archive_filename = filenames[-1]
            else:
                print('WARNING: RII database zip file not found ')
                print('in app data folder')
                print(f'  {appdata_path}')
                print('search in home..')
                filenames = list(Path.home().glob('rii-database-*.zip'))
                if len(filenames) > 0:
                    archive_filename = filenames[-1]
                else:
                    raise Exception('RII database zip file not found.\n'+
                                    'Please download it from url: \n' +
                                    '<https://refractiveindex.info/download/database/>')
        print(f'Using RII database zip file: {archive_filename}')
        self.archive_filename = archive_filename

    def scan(self):
        ''' Scan archive from
            https://refractiveindex.info/download/database/rii-database-2023-10-04.zip
            and output the result in the dict.
            Note: Shelf 'other' is not correctly treated.
        '''
        self.rii_db_items = dict()
        #  0           1               2          3                 4
        # 'database', 'data-nk', 'organic', 'CHBr3 - bromoform', 'Ghosal.yml'
        with zipfile.ZipFile(self.archive_filename) as z:
            for filename in z.namelist():
                if '/' in filename:
                    words = filename.split('/')
                else:
                    words = filename.split('\\')
                if len(words) < 5:
                    continue
                if words[1] == 'data-nk' and len(words[4]) != 0:
                    if words[2] == 'other':
                        key = f'{words[2]}/{words[3]}'
                        if key not in self.rii_db_items:
                            self.rii_db_items[key] = dict()
                        if words[4] not in self.rii_db_items[key]:
                            self.rii_db_items[key][words[4]] = []
                        if words[5][-4:] == '.yml':
                            self.rii_db_items[key][words[4]].append(words[5][:-4])
                    else:
                        if words[2] not in self.rii_db_items:
                            self.rii_db_items[words[2]] = dict()
                        if words[3] not in self.rii_db_items[words[2]]:
                            self.rii_db_items[words[2]][words[3]] = []
                        if words[4][-4:] == '.yml':
                            self.rii_db_items[words[2]][words[3]].append(words[4][:-4])

    def filter_valid(self):
        if self.rii_db_items is None:
            print('No items. Please `_scan()` first')
            return
        filtered_items = dict()
        for shelf in self.rii_db_items:
            for book in self.rii_db_items[shelf]:
                #  print(f'testing {shelf}/{book}')
                for name in self.rii_db_items[shelf][book]:
                    flag = True
                    try:
                        self.select(shelf, book, name)
                    except Exception as e:
                        flag = False
                    if flag:
                        try:
                            self.get_n(300)
                            self.get_k(300)
                            self.get_n(800)
                            self.get_k(800)
                        except Exception as e:
                            flag = False
                    if flag:  # passed
                        if shelf not in filtered_items:
                            filtered_items[shelf] = dict()
                        if book not in filtered_items[shelf]:
                            filtered_items[shelf][book] = []
                        filtered_items[shelf][book].append(name)
        self.rii_db_items = filtered_items

    def print_db_items(self):
        if self.rii_db_items is None:
            print('No items. Please `_scan()` first')
            return
        print('`shelf`')
        print('\t`book`')
        print('\t\t`pages`')
        for shelf in self.rii_db_items:
            print(f'{shelf}')
            for book in self.rii_db_items[shelf]:
                print(f'\t{book}')
                print(f'\t\t{" ".join(self.rii_db_items[shelf][book])}')

    def select(self, shelf, book, name):
        with zipfile.ZipFile(self.archive_filename) as z:
            with z.open(f'database/data-nk/{shelf}/{book}/{name}.yml', 'r') as f:
                data = yaml.safe_load(f)
        i = 0
        type_value = data.get('DATA', [])[i].get('type', '')
        if 'tabulated' not in type_value:
            i = 1
        type_value = data.get('DATA', [])[i].get('type', '')
        if 'tabulated' not in type_value:
            i = 1
        type_value = data.get('DATA', [])[i].get('type', '')
        data_value = data.get('DATA', [])[i].get('data', '')
        lines = data_value.split('\n')
        arrays = [np.array([float(x) for x in s.split()]) for s in lines]

        self.wls = np.array([1000*arr[0] for arr in arrays if arr.size > 0])
        n = np.array([arr[1] for arr in arrays if arr.size > 1])
        k = np.array([arr[2] for arr in arrays if arr.size > 2])
        if len(k) == 0:
            k = np.zeros(len(self.wls))

        self._get_n_interp = interpolate.interp1d(self.wls, n, kind='cubic')
        self._get_k_interp = interpolate.interp1d(self.wls, k, kind='cubic')
        self.__name__ = f'{book}/{name[:5]}'

    def count_books(self):
        i = 0
        for shelf in self.rii_db_items:
            for book in self.rii_db_items[shelf]:
                i += 1
        return i


# ~ def get_file_hash(file_path):
    # ~ hash_md5 = hashlib.md5()
    # ~ with open(file_path, "rb") as f:
        # ~ for chunk in iter(lambda: f.read(4096), b""):
            # ~ hash_md5.update(chunk)
    # ~ return hash_md5.hexdigest()

if __name__ == '__main__':

    riimat = RiiMaterial()

    riimat.scan()
    print(f'Before filtration {riimat.count_books()}')
    riimat.filter_valid()
    print(f'After filtration {riimat.count_books()}')

    riimat.print_db_items()

    riimat.select('main', 'Ag', 'Babar')
    print(riimat)
    riimat.plot()
