import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib
from analiz_sechenia import analiz_sechenia
import os
import subprocess
import glob
import h5py

# Настройки для качественной визуализации
matplotlib.rcParams["font.size"] = 10
matplotlib.rcParams["figure.figsize"] = (12, 8)

# Параметры обнаружения
MIN_DISTANCE = 50

# Параметры окон
WINDOW_SIZE = (128, 128)

# Параметры разрешения из ваших данных
c = 3e8
samp_rate = 720000000.0
range_bandwidth = 600000000.0
wavelength = 0.031228381041666666
azimuth_bandwidth = 240.55034790612274
prf = 522.8468093300061

RANGE_RESOLUTION = c / (2 * range_bandwidth)
AZIMUTH_RESOLUTION = wavelength / (2 * azimuth_bandwidth * wavelength / c)

THRESHOLD_OFFSET_DB = 20

# Папка для результатов
result_folder = "data_20_10_2025"
os.makedirs(result_folder, exist_ok=True)

# Путь к HDF5 файлу
HDF5_FILE_PATH = os.path.join(result_folder, "slc_strip_vv_20250724T113400__F10_24.hdf5")

# ИЗМЕНЕНИЕ 1: Функция для извлечения названия голограммы из пути
def get_hologram_name(file_path):
    """Извлекает название голограммы из пути к файлу"""
    base_name = os.path.basename(file_path)
    name_without_ext = os.path.splitext(base_name)[0]
    return name_without_ext

HOLOGRAM_NAME = get_hologram_name(HDF5_FILE_PATH)

def load_radar_image_from_hdf5(hdf5_path):
    """Загрузка голограммы из HDF5 файла"""
    try:
        with h5py.File(hdf5_path, 'r') as f:
            real_data = None
            imag_data = None
            
            for name in f:
                if 'real' in name.lower():
                    real_data = f[name][:]
                elif 'imag' in name.lower():
                    imag_data = f[name][:]
            
            if real_data is not None and imag_data is not None:
                radar_image_complex = real_data + 1j * imag_data
                radar_image = np.abs(radar_image_complex)
                return radar_image_complex, radar_image, -80
            else:
                raise ValueError("Не найдены real/imag данные в HDF5 файле")
                
    except Exception as e:
        print(f"Ошибка загрузки HDF5: {e}")
        raise

# ... остальные функции без изменений (calculate_noise_threshold, find_targets, extract_target_window, extract_sections, generate_sinc_signal_from_section, convert_to_meters)

def main():
    radar_image_complex, radar_image, noise_power_1 = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    detected_peaks, noise_power_2, SNR_DB_new = find_targets(radar_image, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    
    targets_data = []
    
    for i, target_yx in enumerate(detected_peaks):
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)
        horizontal_section, vertical_section = extract_sections(window)
        
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0]/2)
        
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1]/2)
        
        target_data = {
            'window_linear': window,
            'h_t': t_h, 'h_signal_db': h_signal_db,
            'h_wl': h_results.get('wl'), 'h_wr': h_results.get('wr'),
            'h_width': h_results.get('measured_width', 0),
            'h_pslr': h_results.get('classical_pslr', -80),
            'h_i_pslr': h_results.get('integral_pslr', -80),
            'h_sinc_interp': h_results.get('sinc_interp'),
            'h_t_interp': h_results.get('t_interp'),
            'v_t': t_v, 'v_signal_db': v_signal_db,
            'v_wl': v_results.get('wl'), 'v_wr': v_results.get('wr'),
            'v_width': v_results.get('measured_width', 0),
            'v_pslr': v_results.get('classical_pslr', -80),
            'v_i_pslr': v_results.get('integral_pslr', -80),
            'v_sinc_interp': v_results.get('sinc_interp'),
            'v_t_interp': v_results.get('t_interp'),
        }
        targets_data.append(target_data)

    # Формирование отчета
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Название голограммы: {HOLOGRAM_NAME}\n")
        f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
        f.write(f"Количество целей: {len(detected_peaks)}\n")
        # ... остальные параметры

    # Визуализация (без изменений)
    # ... код визуализации

    # Генерация отчета в Typst
    typ_filename = os.path.join(result_folder, "ozenca.typ")

    # ИЗМЕНЕНИЕ 2: Исправленная нумерация рисунков
    typ_content = f'''#set page(width: auto, height: auto, margin: 1.5cm)
#set text(font: "New Computer Modern", size: 12pt, lang: "ru")
#show heading: set text(weight: "bold")

#align(center)[
#text(size: 24pt, weight: "bold")[Анализ радиолокационного изображения]
]

#align(center)[
#text(size: 12pt)[Голограмма: {HOLOGRAM_NAME}]
]

#align(center)[
#text(size: 12pt)[Размер: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей, количество целей: {len(detected_peaks)}]
]

#figure(
  image("radar_image.png"),
  caption: [Рисунок 1. Радиолокационное изображение с обнаруженными целями]
)'''

    # ИЗМЕНЕНИЕ 3: Правильная нумерация - каждый рисунок получает свой порядковый номер
    figure_counter = 2
    for i, target_data in enumerate(targets_data):
        target_id = i+1
        target_coords = detected_peaks[i]
        target_folder = f"target_{target_id}"
        
        h_width_meters = convert_to_meters(target_data['h_width'], 'horizontal')
        v_width_meters = convert_to_meters(target_data['v_width'], 'vertical')
        
        # Первый grid - окна цели (2 рисунка)
        typ_content += f'''
#pagebreak()
#align(center)[
#text(size: 18pt, weight: "bold")[Цель №{target_id}]
]

#align(center)[
#text(size: 12pt)[Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}]
]

// Визуализация окна цели
#figure(
  grid(
    columns: 2,
    gutter: 2cm,
    [
      #figure(
        image("{target_folder}/target_{target_id}_linear.png"),
        caption: [Рисунок {figure_counter}. Окно цели - линейный масштаб]
      )
    ],
    [
      #figure(
        image("{target_folder}/target_{target_id}_db.png"),
        caption: [Рисунок {figure_counter+1}. Окно цели - логарифмический масштаб]
      )
    ]
  )
)'''
        figure_counter += 2  # Увеличиваем на 2 для двух рисунков

        # Второй grid - сечения (2 рисунка)
        typ_content += f'''
// Сечения цели
#figure(
  grid(
    columns: 2,
    gutter: 2cm,
    [
      #figure(
        image("{target_folder}/target_{target_id}_horizontal.png"),
        caption: [Рисунок {figure_counter}. Сечение по дальности]
      )
    ],
    [
      #figure(
        image("{target_folder}/target_{target_id}_vertical.png"),
        caption: [Рисунок {figure_counter+1}. Сечение по азимуту]
      )
    ]
  )
)'''
        figure_counter += 2  # Увеличиваем на 2 для двух рисунков

        # Таблицы с параметрами
        typ_content += f'''
#figure(
  grid(
    columns: 2,
    gutter: 3cm,
    [
      #table(
        columns: 2,
        align: center,
        stroke: (x: 0.5pt, y: 0.5pt),
        inset: 5pt,
        [*Параметр сечения по дальности*], [*Значение*],
        [Ширина главного лепестка], [{target_data['h_width']:.4f} отсч. ({h_width_meters:.4f} м)],
        [Максимальный УБЛ], [{target_data['h_pslr']:.2f} дБ],
        [Интегральный УБЛ], [{target_data['h_i_pslr']:.2f} дБ],
      )
    ],
    [
      #table(
        columns: 2,
        align: center,
        stroke: (x: 0.5pt, y: 0.5pt),
        inset: 5pt,
        [*Параметр сечения по азимуту*], [*Значение*],
        [Ширина главного лепестка], [{target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)],
        [Максимальный УБЛ], [{target_data['v_pslr']:.2f} дБ],
        [Интегральный УБЛ], [{target_data['v_i_pslr']:.2f} дБ],
      )
    ]
  )
)'''
        # Для таблиц не увеличиваем figure_counter, так как они не считаются рисунками

    with open(typ_filename, 'w', encoding='utf-8') as typ_file:
        typ_file.write(typ_content)

    # Компиляция PDF и очистка (без изменений)
    # ... код компиляции

if __name__ == "__main__":
    main()
