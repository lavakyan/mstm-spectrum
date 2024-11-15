import yaml
import numpy as np
from mstm_studio.mstm_spectrum import Material

def import_from_file_plot(filename):
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)

    # Извлечение данных из YAML
    # Разделение строки данных на массивы для удобства
    data_points = data['DATA'][0]['data'].strip().splitlines()
    wls = []
    n_values = []
    k_values = []

    for line in data_points:
        wl, n, k = map(float, line.split())
        # Преобразуем длину волны из микрон в нанометры
        wls.append(wl * 1000)  
        n_values.append(n)
        k_values.append(k)

    # Преобразование списков в numpy массивы
    wls = np.array(wls)
    n_values = np.array(n_values)
    k_values = np.array(k_values)
    # Разделяем строку по символу '\'
    parts = filename.split('\\')

    # Получаем нужные части строки
    element = parts[-2]
    compound = parts[-1].split('.')[0] 

    # Формируем итоговую строку
    name = f"{element}({compound})"
    # Создание экземпляра класса Material
    return Material(
        file_name= name,  # Название материала
        wls=wls,
        nk=n_values + 1j * k_values  # Комплексное значение nk из n и k
    )