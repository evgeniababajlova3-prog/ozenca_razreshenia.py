import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib
from analiz_sechenia import analiz_sechenia
import os
import subprocess
import glob
import h5py  # ИЗМЕНЕНИЕ 1: Добавлен импорт h5py для работы с HDF5

# Настройки для качественной визуализации
matplotlib.rcParams["font.size"] = 10
matplotlib.rcParams["figure.figsize"] = (12, 8)

# Параметры обнаружения
MIN_DISTANCE = 50

# Параметры окон
WINDOW_SIZE = (128, 128)

# ИЗМЕНЕНИЕ 2: Оптимальные параметры разрешения из ваших данных съёмки
# Параметры из JSON файла
c = 3e8  # скорость света, м/с

# По дальности (из ваших параметров)
samp_rate = 720000000.0  # Гц
range_bandwidth = 600000000.0  # Гц
RANGE_RESOLUTION = c / (2 * range_bandwidth)  # метров на отсчет

# По азимуту (из ваших параметров)  
wavelength = 0.031228381041666666  # м
azimuth_bandwidth = 240.55034790612274  # Гц
prf = 522.8468093300061  # Гц

# Оптимальная формула для азимутального разрешения в SAR
# Разрешение = скорость света / (2 * доплеровская полоса * наклонная дальность / скорость)
# Упрощенный вариант: используем стандартную формулу через доплеровскую полосу
AZIMUTH_RESOLUTION = wavelength / (2 * azimuth_bandwidth * wavelength / c)  # метров на отсчет
# Альтернативная упрощенная формула (часто используемая):
# AZIMUTH_RESOLUTION = prf / azimuth_bandwidth * (c / (2 * samp_rate))

THRESHOLD_OFFSET_DB = 20

# Папка для результатов
result_folder = "data_20_10_2025"
os.makedirs(result_folder, exist_ok=True)

# Путь к HDF5 файлу ВНУТРИ папки с результатами
HDF5_FILE_PATH = os.path.join(result_folder, "slc_strip_vv_20250724T113400__F10_24.hdf5")

# ИЗМЕНЕНИЕ 3: Добавлена функция загрузки HDF5 вместо генерации РЛИ
def load_radar_image_from_hdf5(hdf5_path):
    """
    Загрузка голограммы из HDF5 файла (замена generate_radar_image)
    """
    try:
        with h5py.File(hdf5_path, 'r') as f:
            # Ищем real и imag части
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

############################################
##
# МОДУЛЬ 2: ОБНАРУЖЕНИЕ ЦЕЛЕЙ (без изменений)
############################################
##

def calculate_noise_threshold(radar_image_complex, offset_db):
    """Расчет порога на основе мощности шума в области с распределением Релея и равномерными фазами."""
    noise_window, amplitudes, phases, match_quality = find_rayleigh_uniform_region(radar_image_complex)
    noise_power = np.mean(amplitudes ** 2)
    noise_power_db = 10 * np.log10(noise_power)
    noise_rms = np.sqrt(noise_power)
    threshold_db = 20 * np.log10(noise_rms) + offset_db
    threshold_linear = 10 ** (threshold_db / 20)
    return threshold_linear, noise_power_db

def find_targets(radar_image, min_distance, threshold_offset_db=10):
    """Обнаружение целей."""
    threshold, noise_power_2 = calculate_noise_threshold(radar_image, threshold_offset_db)
    local_max = ndimage.maximum_filter(radar_image, size=min_distance) == radar_image
    above_threshold = radar_image > threshold
    detected = local_max & above_threshold
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))
    SNR_DB_new = 10 * np.log10(np.max(np.abs(radar_image))) - noise_power_2
    return peaks_coords, noise_power_2, SNR_DB_new

############################################
##
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (с минимальными изменениями)
############################################
##

def extract_target_window(radar_image, center_yx, window_size):
    """Выделяет окно вокруг цели."""
    center_y, center_x = center_yx
    win_h, win_w = window_size
    half_h, half_w = win_h // 2, win_w // 2

    y_start = max(0, center_y - half_h)
    y_end = min(radar_image.shape[0], center_y + half_h)
    x_start = max(0, center_x - half_w)
    x_end = min(radar_image.shape[1], center_x + half_w)

    window = radar_image[y_start:y_end, x_start:x_end]
    return window

def extract_sections(window):
    """Извлекает горизонтальное и вертикальное сечения из окна."""
    center_y, center_x = window.shape[0] // 2, window.shape[1] // 2
    horizontal_section = window[center_y, :]  # Центральная строка
    vertical_section = window[:, center_x]    # Центральный столбец
    return horizontal_section, vertical_section

def generate_sinc_signal_from_section(section, window_size):
    """Преобразует сечение в формат для анализа sinc-функции."""
    section_linear = np.abs(section)
    section_norm = section_linear / np.max(section_linear)
    section_db = 20 * np.log10(section_norm)
    
    # ИЗМЕНЕНИЕ 4: ось от 0 до 128 вместо -64 до +64
    t = np.linspace(0, window_size, len(section))
    
    return t, section_db, section_norm

# ИЗМЕНЕНИЕ 5: Добавлена функция конвертации в метры
def convert_to_meters(width_samples, direction):
    """Конвертирует ширину из отсчетов в метры."""
    if direction == 'horizontal':
        return width_samples * RANGE_RESOLUTION
    elif direction == 'vertical':
        return width_samples * AZIMUTH_RESOLUTION
    else:
        return width_samples

############################################
##
# ОСНОВНАЯ ФУНКЦИЯ (с изменениями визуализации)
############################################
##

def main():
    # ИЗМЕНЕНИЕ 6: Замена генерации на загрузку из HDF5
    radar_image_complex, radar_image, noise_power_1 = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    
    # Обнаружение целей (без изменений)
    detected_peaks, noise_power_2, SNR_DB_new = find_targets(radar_image, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    
    targets_data = []
    
    for i, target_yx in enumerate(detected_peaks):
        # Выделение окна
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)
        
        # Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window)
        
        # Анализ горизонтального сечения
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0]/2)
        
        # Анализ вертикального сечения
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1]/2)
        
        # Формируем данные для отчета
        target_data = {
            'window_linear': window,
            'h_t': t_h,
            'h_signal_db': h_signal_db,
            'h_wl': h_results.get('wl'),
            'h_wr': h_results.get('wr'),
            'h_width': h_results.get('measured_width', 0),
            'h_pslr': h_results.get('classical_pslr', -80),
            'h_i_pslr': h_results.get('integral_pslr', -80),
            'h_sinc_interp': h_results.get('sinc_interp'),
            'h_t_interp': h_results.get('t_interp'),
            'v_t': t_v,
            'v_signal_db': v_signal_db,
            'v_wl': v_results.get('wl'),
            'v_wr': v_results.get('wr'),
            'v_width': v_results.get('measured_width', 0),
            'v_pslr': v_results.get('classical_pslr', -80),
            'v_i_pslr': v_results.get('integral_pslr', -80),
            'v_sinc_interp': v_results.get('sinc_interp'),
            'v_t_interp': v_results.get('t_interp'),
        }
        
        targets_data.append(target_data)

    # Формирование отчета
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
        f.write(f"Количество целей: {len(detected_peaks)}\n")
        f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
        f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
        f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")
        # ИЗМЕНЕНИЕ 7: Добавлена детальная информация о параметрах съёмки
        f.write(f"\nПараметры съёмки:\n")
        f.write(f"Частота дискретизации: {samp_rate/1e6:.1f} МГц\n")
        f.write(f"Несущая частота: {9.6e9/1e9:.1f} ГГц\n")
        f.write(f"Длина волны: {wavelength:.6f} м\n")
        f.write(f"Полоса по дальности: {range_bandwidth/1e6:.1f} МГц\n")
        f.write(f"Полоса по азимуту: {azimuth_bandwidth:.1f} Гц\n")
        f.write(f"PRF: {prf:.1f} Гц\n")
        f.write(f"Разрешение по дальности: {RANGE_RESOLUTION:.6f} м/отсчет\n")
        f.write(f"Разрешение по азимуту: {AZIMUTH_RESOLUTION:.6f} м/отсчет\n")

    # Сохраняем основное РЛИ с целями
    plt.figure(figsize=(9, 8))
    plt.imshow(radar_image, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Амплитуда')
    plt.title('Радиолокационное изображение с обнаруженными целями')
    plt.xlabel('Координаты по дальности')
    plt.ylabel('Координаты по азимуту')

    # Отмечаем обнаруженные цели
    for i, (y, x) in enumerate(detected_peaks):
        plt.plot(x, y, 's', markersize=14, markeredgewidth=2,
                markeredgecolor='red', markerfacecolor='none', linestyle='none',
                label='Цели' if i == 0 else "")
    
    plt.legend()
    radar_image_path = os.path.join(result_folder, "radar_image.png")
    plt.savefig(radar_image_path, dpi=150, bbox_inches='tight')
    plt.close()

    # Сохраняем данные и графики для каждой цели
    for i, target_data in enumerate(targets_data):
        target_id = i+1
        target_coords = detected_peaks[i]

        # Создаем папку для цели
        target_folder = os.path.join(result_folder, f"target_{target_id}")
        os.makedirs(target_folder, exist_ok=True)

        # ИЗМЕНЕНИЕ 8: Рассчитываем ширину в метрах
        h_width_meters = convert_to_meters(target_data['h_width'], 'horizontal')
        v_width_meters = convert_to_meters(target_data['v_width'], 'vertical')

        # Сохраняем параметры цели
        with open(os.path.join(target_folder, f"target_{target_id}_params.txt"), 'w', encoding='utf-8') as f:
            f.write(f"Цель №{target_id}\n")
            f.write(f"Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}\n\n")

            f.write("Сечение по дальности:\n")
            # ИЗМЕНЕНИЕ 9: Добавлены метры рядом с отсчетами
            f.write(f"Ширина главного лепестка: {target_data['h_width']:.4f} отсч. ({h_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['h_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['h_i_pslr']:.2f} дБ\n\n")

            f.write("Сечение по азимуту:\n")
            # ИЗМЕНЕНИЕ 10: Добавлены метры рядом с отсчетами
            f.write(f"Ширина главного лепестка: {target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['v_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['v_i_pslr']:.2f} дБ\n")

        # Сохраняем массивы данных
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_window_linear.txt"),
                  target_data['window_linear'], fmt='%6f')
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_horizontal_section.txt"),
                  np.column_stack([target_data['h_t'], target_data['h_signal_db']]), fmt='%6f')
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_vertical_section.txt"),
                  np.column_stack([target_data['v_t'], target_data['v_signal_db']]), fmt='%6f')

        # Окно в линейном масштабе
        plt.figure(figsize=(6, 4))
        plt.imshow(target_data['window_linear'], cmap='hot', interpolation='nearest')
        # ИЗМЕНЕНИЕ 11: Исправлена подпись - убрано "(дБ)" для линейного масштаба
        plt.colorbar(label='Амплитуда')
        plt.title(f'Цель {target_id} - Линейный масштаб')
        plt.xlabel('Дальность')
        plt.ylabel('Азимут')
        linear_path = os.path.join(target_folder, f"target_{target_id}_linear.png")
        plt.savefig(linear_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Окно в логарифмическом масштабе
        plt.figure(figsize=(6, 4))
        window_db = 20 * np.log10(np.abs(target_data['window_linear']))
        plt.imshow(window_db, cmap='hot', interpolation='nearest')
        plt.colorbar(label='Амплитуда (дБ)')
        plt.title(f'Цель {target_id} - Логарифмический масштаб')
        plt.xlabel('Дальность')
        plt.ylabel('Азимут')
        db_path = os.path.join(target_folder, f"target_{target_id}_db.png")
        plt.savefig(db_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Горизонтальное сечение
        plt.figure(figsize=(14, 7))
        
        # Исходный и интерполированный сигнал
        plt.plot(target_data['h_t'], target_data['h_signal_db'], 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
        plt.plot(target_data['h_t_interp'], target_data['h_sinc_interp'], 'g-', linewidth=1.5, label='Интерполированный сигнал', alpha=0.8)

        # ИЗМЕНЕНИЕ 12: Показываем ширину главного лепестка (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['h_wl'] is not None and target_data['h_wr'] is not None:
            plt.axvline(x=target_data['h_wl'], color='red', linestyle='--', alpha=0.7)
            plt.axvline(x=target_data['h_wr'], color='red', linestyle='--', alpha=0.7)

        # ИЗМЕНЕНИЕ 13: Показываем максимальный УБЛ (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['h_pslr'] > -80:
            pslr_level = -target_data['h_pslr']
            plt.axhline(y=pslr_level, color='purple', linestyle='-', alpha=0.7)

        plt.xlabel('Дальность, отсчеты')
        plt.ylabel('Амплитуда (дБ)')
        plt.title(f'Анализ сечения по дальности - цель {target_id}')
        plt.legend(fontsize=16)
        plt.grid(True, alpha=0.5, linewidth=0.8)
        plt.ylim(-50, 5)
        plt.tight_layout()
        h_section_path = os.path.join(target_folder, f"target_{target_id}_horizontal.png")
        plt.savefig(h_section_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Вертикальное сечение
        plt.figure(figsize=(14, 7))
        
        # Исходный и интерполированный сигнал
        plt.plot(target_data['v_t'], target_data['v_signal_db'], 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
        plt.plot(target_data['v_t_interp'], target_data['v_sinc_interp'], 'g-', linewidth=1.5, label='Интерполированный сигнал', alpha=0.8)

        # ИЗМЕНЕНИЕ 14: Показываем ширину главного лепестка (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['v_wl'] is not None and target_data['v_wr'] is not None:
            plt.axvline(x=target_data['v_wl'], color='red', linestyle='--', alpha=0.7)
            plt.axvline(x=target_data['v_wr'], color='red', linestyle='--', alpha=0.7)

        # ИЗМЕНЕНИЕ 15: Показываем максимальный УБЛ (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['v_pslr'] > -80:
            pslr_level = -target_data['v_pslr']
            plt.axhline(y=pslr_level, color='purple', linestyle='-', alpha=0.7)

        plt.xlabel('Азимут, отсчеты')
        plt.ylabel('Амплитуда (дБ)')
        plt.title(f'Анализ сечения по азимуту - цель {target_id}')
        plt.legend(fontsize=16)
        plt.grid(True, alpha=0.5, linewidth=0.8)
        plt.ylim(-50, 5)
        plt.tight_layout()
        v_section_path = os.path.join(target_folder, f"target_{target_id}_vertical.png")
        plt.savefig(v_section_path, dpi=150, bbox_inches='tight')
        plt.close()

    # Генерация отчета в Typst
    typ_filename = os.path.join(result_folder, "ozenca.typ")

    # ИЗМЕНЕНИЕ 16: Генерируем typst-код с правильной нумерацией (последовательной)
    typ_content = f'''#set page(width: auto, height: auto, margin: 1.5cm)
#set text(font: "New Computer Modern", size: 12pt, lang: "ru")
#show heading: set text(weight: "bold")

#align(center)[
#text(size: 24pt, weight: "bold")[Анализ радиолокационного изображения]
]

#align(center)[
#text(size: 12pt)[Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей, количество целей: {len(detected_peaks)}]
]

#figure(
  image("radar_image.png"),
  caption: [Рисунок 1. Радиолокационное изображение с обнаруженными целями]
)'''

    # ИЗМЕНЕНИЕ 17: Добавляем данные по каждой цели с последовательной нумерацией
    figure_counter = 2
    for i, target_data in enumerate(targets_data):
        target_id = i+1
        target_coords = detected_peaks[i]
        target_folder = f"target_{target_id}"
        
        # ИЗМЕНЕНИЕ 18: Рассчитываем ширину в метрах для таблицы
        h_width_meters = convert_to_meters(target_data['h_width'], 'horizontal')
        v_width_meters = convert_to_meters(target_data['v_width'], 'vertical')
        
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
)

// Сечения цели
#figure(
  grid(
    columns: 2,
    gutter: 2cm,
    [
      #figure(
        image("{target_folder}/target_{target_id}_horizontal.png"),
        caption: [Рисунок {figure_counter+2}. Сечение по дальности]
      )
    ],
    [
      #figure(
        image("{target_folder}/target_{target_id}_vertical.png"),
        caption: [Рисунок {figure_counter+3}. Сечение по азимуту]
      )
    ]
  )
)

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
        // ИЗМЕНЕНИЕ 19: Добавлены метры в таблицу
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
        // ИЗМЕНЕНИЕ 20: Добавлены метры в таблицу
        [Ширина главного лепестка], [{target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)],
        [Максимальный УБЛ], [{target_data['v_pslr']:.2f} дБ],
        [Интегральный УБЛ], [{target_data['v_i_pslr']:.2f} дБ],
      )
    ]
  )
)'''
        figure_counter += 4  # ИЗМЕНЕНИЕ 21: Увеличиваем счетчик на 4 для каждой цели

    with open(typ_filename, 'w', encoding='utf-8') as typ_file:
        typ_file.write(typ_content)

    # Компиляция PDF
    typ_path = os.path.join(result_folder, "ozenca.typ")
    pdf_path = typ_path.replace('.typ', '.pdf')
    typst_path = r'C:\Users\Gisich_AV\Desktop\typst\typst\typst.exe'
    
    if os.path.exists(typst_path):
        # ИЗМЕНЕНИЕ 22: Убраны выводы в консоль
        result = subprocess.run([typst_path, 'compile', typ_path, pdf_path], 
                              capture_output=True, text=True)

    # Очистка временных файлов
    pdf_files = glob.glob(os.path.join(result_folder, "*.pdf"))
    for file in glob.glob(os.path.join(result
# ИЗМЕНЕНИЕ 3: Добавлена функция загрузки HDF5 вместо генерации РЛИ
def load_radar_image_from_hdf5(hdf5_path):
    """
    Загрузка голограммы из HDF5 файла (замена generate_radar_image)
    """
    try:
        with h5py.File(hdf5_path, 'r') as f:
            # Ищем real и imag части
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

############################################
##
# МОДУЛЬ 2: ОБНАРУЖЕНИЕ ЦЕЛЕЙ (без изменений)
############################################
##

def calculate_noise_threshold(radar_image_complex, offset_db):
    """Расчет порога на основе мощности шума в области с распределением Релея и равномерными фазами."""
    noise_window, amplitudes, phases, match_quality = find_rayleigh_uniform_region(radar_image_complex)
    noise_power = np.mean(amplitudes ** 2)
    noise_power_db = 10 * np.log10(noise_power)
    noise_rms = np.sqrt(noise_power)
    threshold_db = 20 * np.log10(noise_rms) + offset_db
    threshold_linear = 10 ** (threshold_db / 20)
    return threshold_linear, noise_power_db

def find_targets(radar_image, min_distance, threshold_offset_db=10):
    """Обнаружение целей."""
    threshold, noise_power_2 = calculate_noise_threshold(radar_image, threshold_offset_db)
    local_max = ndimage.maximum_filter(radar_image, size=min_distance) == radar_image
    above_threshold = radar_image > threshold
    detected = local_max & above_threshold
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))
    SNR_DB_new = 10 * np.log10(np.max(np.abs(radar_image))) - noise_power_2
    return peaks_coords, noise_power_2, SNR_DB_new

############################################
##
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (с минимальными изменениями)
############################################
##

def extract_target_window(radar_image, center_yx, window_size):
    """Выделяет окно вокруг цели."""
    center_y, center_x = center_yx
    win_h, win_w = window_size
    half_h, half_w = win_h // 2, win_w // 2

    y_start = max(0, center_y - half_h)
    y_end = min(radar_image.shape[0], center_y + half_h)
    x_start = max(0, center_x - half_w)
    x_end = min(radar_image.shape[1], center_x + half_w)

    window = radar_image[y_start:y_end, x_start:x_end]
    return window

def extract_sections(window):
    """Извлекает горизонтальное и вертикальное сечения из окна."""
    center_y, center_x = window.shape[0] // 2, window.shape[1] // 2
    horizontal_section = window[center_y, :]  # Центральная строка
    vertical_section = window[:, center_x]    # Центральный столбец
    return horizontal_section, vertical_section

def generate_sinc_signal_from_section(section, window_size):
    """Преобразует сечение в формат для анализа sinc-функции."""
    section_linear = np.abs(section)
    section_norm = section_linear / np.max(section_linear)
    section_db = 20 * np.log10(section_norm)
    
    # ИЗМЕНЕНИЕ 4: ось от 0 до 128 вместо -64 до +64
    t = np.linspace(0, window_size, len(section))
    
    return t, section_db, section_norm

# ИЗМЕНЕНИЕ 5: Добавлена функция конвертации в метры
def convert_to_meters(width_samples, direction):
    """Конвертирует ширину из отсчетов в метры."""
    if direction == 'horizontal':
        return width_samples * RANGE_RESOLUTION
    elif direction == 'vertical':
        return width_samples * AZIMUTH_RESOLUTION
    else:
        return width_samples

############################################
##
# ОСНОВНАЯ ФУНКЦИЯ (с изменениями визуализации)
############################################
##

def main():
    # ИЗМЕНЕНИЕ 6: Замена генерации на загрузку из HDF5
    radar_image_complex, radar_image, noise_power_1 = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    
    # Обнаружение целей (без изменений)
    detected_peaks, noise_power_2, SNR_DB_new = find_targets(radar_image, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    
    targets_data = []
    
    for i, target_yx in enumerate(detected_peaks):
        # Выделение окна
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)
        
        # Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window)
        
        # Анализ горизонтального сечения
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0]/2)
        
        # Анализ вертикального сечения
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1]/2)
        
        # Формируем данные для отчета
        target_data = {
            'window_linear': window,
            'h_t': t_h,
            'h_signal_db': h_signal_db,
            'h_wl': h_results.get('wl'),
            'h_wr': h_results.get('wr'),
            'h_width': h_results.get('measured_width', 0),
            'h_pslr': h_results.get('classical_pslr', -80),
            'h_i_pslr': h_results.get('integral_pslr', -80),
            'h_sinc_interp': h_results.get('sinc_interp'),
            'h_t_interp': h_results.get('t_interp'),
            'v_t': t_v,
            'v_signal_db': v_signal_db,
            'v_wl': v_results.get('wl'),
            'v_wr': v_results.get('wr'),
            'v_width': v_results.get('measured_width', 0),
            'v_pslr': v_results.get('classical_pslr', -80),
            'v_i_pslr': v_results.get('integral_pslr', -80),
            'v_sinc_interp': v_results.get('sinc_interp'),
            'v_t_interp': v_results.get('t_interp'),
        }
        
        targets_data.append(target_data)

    # Формирование отчета
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
        f.write(f"Количество целей: {len(detected_peaks)}\n")
        f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
        f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
        f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")
        # ИЗМЕНЕНИЕ 7: Добавлена детальная информация о параметрах съёмки
        f.write(f"\nПараметры съёмки:\n")
        f.write(f"Частота дискретизации: {samp_rate/1e6:.1f} МГц\n")
        f.write(f"Несущая частота: {9.6e9/1e9:.1f} ГГц\n")
        f.write(f"Длина волны: {wavelength:.6f} м\n")
        f.write(f"Полоса по дальности: {range_bandwidth/1e6:.1f} МГц\n")
        f.write(f"Полоса по азимуту: {azimuth_bandwidth:.1f} Гц\n")
        f.write(f"PRF: {prf:.1f} Гц\n")
        f.write(f"Разрешение по дальности: {RANGE_RESOLUTION:.6f} м/отсчет\n")
        f.write(f"Разрешение по азимуту: {AZIMUTH_RESOLUTION:.6f} м/отсчет\n")

    # Сохраняем основное РЛИ с целями
    plt.figure(figsize=(9, 8))
    plt.imshow(radar_image, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Амплитуда')
    plt.title('Радиолокационное изображение с обнаруженными целями')
    plt.xlabel('Координаты по дальности')
    plt.ylabel('Координаты по азимуту')

    # Отмечаем обнаруженные цели
    for i, (y, x) in enumerate(detected_peaks):
        plt.plot(x, y, 's', markersize=14, markeredgewidth=2,
                markeredgecolor='red', markerfacecolor='none', linestyle='none',
                label='Цели' if i == 0 else "")
    
    plt.legend()
    radar_image_path = os.path.join(result_folder, "radar_image.png")
    plt.savefig(radar_image_path, dpi=150, bbox_inches='tight')
    plt.close()

    # Сохраняем данные и графики для каждой цели
    for i, target_data in enumerate(targets_data):
        target_id = i+1
        target_coords = detected_peaks[i]

        # Создаем папку для цели
        target_folder = os.path.join(result_folder, f"target_{target_id}")
        os.makedirs(target_folder, exist_ok=True)

        # ИЗМЕНЕНИЕ 8: Рассчитываем ширину в метрах
        h_width_meters = convert_to_meters(target_data['h_width'], 'horizontal')
        v_width_meters = convert_to_meters(target_data['v_width'], 'vertical')

        # Сохраняем параметры цели
        with open(os.path.join(target_folder, f"target_{target_id}_params.txt"), 'w', encoding='utf-8') as f:
            f.write(f"Цель №{target_id}\n")
            f.write(f"Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}\n\n")

            f.write("Сечение по дальности:\n")
            # ИЗМЕНЕНИЕ 9: Добавлены метры рядом с отсчетами
            f.write(f"Ширина главного лепестка: {target_data['h_width']:.4f} отсч. ({h_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['h_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['h_i_pslr']:.2f} дБ\n\n")

            f.write("Сечение по азимуту:\n")
            # ИЗМЕНЕНИЕ 10: Добавлены метры рядом с отсчетами
            f.write(f"Ширина главного лепестка: {target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['v_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['v_i_pslr']:.2f} дБ\n")

        # Сохраняем массивы данных
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_window_linear.txt"),
                  target_data['window_linear'], fmt='%6f')
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_horizontal_section.txt"),
                  np.column_stack([target_data['h_t'], target_data['h_signal_db']]), fmt='%6f')
        np.savetxt(os.path.join(target_folder, f"target_{target_id}_vertical_section.txt"),
                  np.column_stack([target_data['v_t'], target_data['v_signal_db']]), fmt='%6f')

        # Окно в линейном масштабе
        plt.figure(figsize=(6, 4))
        plt.imshow(target_data['window_linear'], cmap='hot', interpolation='nearest')
        # ИЗМЕНЕНИЕ 11: Исправлена подпись - убрано "(дБ)" для линейного масштаба
        plt.colorbar(label='Амплитуда')
        plt.title(f'Цель {target_id} - Линейный масштаб')
        plt.xlabel('Дальность')
        plt.ylabel('Азимут')
        linear_path = os.path.join(target_folder, f"target_{target_id}_linear.png")
        plt.savefig(linear_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Окно в логарифмическом масштабе
        plt.figure(figsize=(6, 4))
        window_db = 20 * np.log10(np.abs(target_data['window_linear']))
        plt.imshow(window_db, cmap='hot', interpolation='nearest')
        plt.colorbar(label='Амплитуда (дБ)')
        plt.title(f'Цель {target_id} - Логарифмический масштаб')
        plt.xlabel('Дальность')
        plt.ylabel('Азимут')
        db_path = os.path.join(target_folder, f"target_{target_id}_db.png")
        plt.savefig(db_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Горизонтальное сечение
        plt.figure(figsize=(14, 7))
        
        # Исходный и интерполированный сигнал
        plt.plot(target_data['h_t'], target_data['h_signal_db'], 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
        plt.plot(target_data['h_t_interp'], target_data['h_sinc_interp'], 'g-', linewidth=1.5, label='Интерполированный сигнал', alpha=0.8)

        # ИЗМЕНЕНИЕ 12: Показываем ширину главного лепестка (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['h_wl'] is not None and target_data['h_wr'] is not None:
            plt.axvline(x=target_data['h_wl'], color='red', linestyle='--', alpha=0.7)
            plt.axvline(x=target_data['h_wr'], color='red', linestyle='--', alpha=0.7)

        # ИЗМЕНЕНИЕ 13: Показываем максимальный УБЛ (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['h_pslr'] > -80:
            pslr_level = -target_data['h_pslr']
            plt.axhline(y=pslr_level, color='purple', linestyle='-', alpha=0.7)

        plt.xlabel('Дальность, отсчеты')
        plt.ylabel('Амплитуда (дБ)')
        plt.title(f'Анализ сечения по дальности - цель {target_id}')
        plt.legend(fontsize=16)
        plt.grid(True, alpha=0.5, linewidth=0.8)
        plt.ylim(-50, 5)
        plt.tight_layout()
        h_section_path = os.path.join(target_folder, f"target_{target_id}_horizontal.png")
        plt.savefig(h_section_path, dpi=150, bbox_inches='tight')
        plt.close()

        # Вертикальное сечение
        plt.figure(figsize=(14, 7))
        
        # Исходный и интерполированный сигнал
        plt.plot(target_data['v_t'], target_data['v_signal_db'], 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
        plt.plot(target_data['v_t_interp'], target_data['v_sinc_interp'], 'g-', linewidth=1.5, label='Интерполированный сигнал', alpha=0.8)

        # ИЗМЕНЕНИЕ 14: Показываем ширину главного лепестка (БЕЗ ЛЕГЕНДЫ - убраны label)
        if target_data['v_wl'] is not None and target_data['v_wr'] is not None:
            plt.axvline(x=target_data['v_wl'], color='red', linestyle='--', alpha=0.7)
            plt.axvline(x=target_data['v_wr'], color
