import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib
from analiz_sechenia import analiz_sechenia
import os
import subprocess

# Настройки для качественной визуализации
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['figure.figsize'] = (12, 8)

# Априорные параметры РЛИ
IMAGE_SIZE = (2048, 2048)  # Размер радиолокационного изображения
TARGETS = [(1500, 100), (1350, 700)]  # Координаты целей (x, y)
SINC_WIDTH_X = 2.0  # Ширина sinc-функции в горизонтальном сечении
SINC_WIDTH_Y = 3.0 # Ширина sinc-функции в вертикальном сечении

# Параметры обнаружения
MIN_DISTANCE = 50

# Параметры окон
WINDOW_SIZE = (128, 128)  # Размер окна анализа вокруг цели
SINC_WIDTH = 1.0
DISCR_PARAM = 5
TIME_RANGE = 10.0

RANGE_RESOLUTION = 1.0  # метров на отсчет по дальности
AZIMUTH_RESOLUTION = 1.0  # метров на отсчет по азимуту

result_folder = "data_20_10_2025"
os.makedirs(result_folder, exist_ok=True)

###############################################################################
# МОДУЛЬ 1: ГЕНЕРАЦИЯ РАДИОЛОКАЦИОННОГО ИЗОБРАЖЕНИЯ
###############################################################################

def generate_2d_sinc(x0, y0, size, lobe_width_X, lobe_width_Y, discr_param):
    """
    Генерирует двумерную sinc-функцию как произведение двух одномерных sinc.
    """
    # Создаем координатные сетки
    y, x = np.ogrid[0:size[0], 0:size[1]]

    # Вычисляем нормализованные расстояния от центра
    x_norm = (x - x0) / (lobe_width_X * discr_param)
    y_norm = (y - y0) / (lobe_width_Y * discr_param)

    # Создаем двумерную sinc как произведение двух одномерных
    sinc_x = np.sinc(x_norm)
    sinc_y = np.sinc(y_norm)
    sinc_2d = sinc_x * sinc_y

    return sinc_2d


def generate_radar_image(targets, image_size, lobe_width_X, lobe_width_Y, discr_param, noise_level_db):
    """Генерирует радиолокационное изображение с реалистичным шумом."""
    radar_image = np.zeros(image_size)

    for target in targets:
        x, y = target
        sinc_target = generate_2d_sinc(x, y, image_size, lobe_width_X, lobe_width_Y, discr_param)
        radar_image += sinc_target

    # Добавляем реалистичный шум
    peak_signal = np.max(radar_image)
    noise_power = peak_signal * 10 ** (noise_level_db / 20)
    noise_real = np.random.normal(0, noise_power, image_size)
    noise_imag = np.random.normal(0, noise_power, image_size)
    complex_noise = noise_real + 1j * noise_imag

    # Преобразуем в комплексный сигнал с шумом
    radar_image_complex = radar_image + complex_noise
    radar_image = np.abs(radar_image_complex)

    return radar_image


###############################################################################
# МОДУЛЬ 2: ОБНАРУЖЕНИЕ ЦЕЛЕЙ
###############################################################################

def calculate_noise_threshold(radar_image, window_size=(32, 32), x_db=10):
    """
    расчет порога по шумовой области.

    Args:
        radar_image: радиолокационное изображение
        window_size: размер окна для анализа шума
        x_db: добавка в дБ к медианному уровню шума
    """
    h, w = radar_image.shape
    window_h, window_w = window_size

    # 1. Ищем углы изображения как заведомо шумовые области
    corners = [
        radar_image[0:window_h, 0:window_w],  # левый верхний
        radar_image[0:window_h, w - window_w:w],  # правый верхний
        radar_image[h - window_h:h, 0:window_w],  # левый нижний
        radar_image[h - window_h:h, w - window_w:w]  # правый нижний
    ]

    # 2. Добавляем несколько случайных окон, но с проверкой на "нецелевые" характеристики
    additional_windows = []
    for i in range(20):
        y_start = np.random.randint(0, h - window_h)
        x_start = np.random.randint(0, w - window_w)

        window = radar_image[y_start:y_start + window_h, x_start:x_start + window_w]

        # Критерии для шумового окна:
        # - низкая максимальная амплитуда (меньше 10% от максимума изображения)
        # - низкое отношение максимума к среднему
        max_val = np.max(window)
        mean_val = np.mean(window)

        if (max_val < 0.1 * np.max(radar_image) and
                max_val / mean_val < 3.0):  # в равномерном шуме это отношение небольшое
            additional_windows.append(window)

    # 3. Собираем все кандидаты в шумовые области
    all_noise_windows = corners + additional_windows


        # Вычисляем медиану средних значений по окнам (медиана устойчива к выбросам)
    window_means = [np.mean(window) for window in all_noise_windows]
    noise_estimate = np.median(window_means)

    # 4. Переводим в дБ и добавляем X дБ
    noise_estimate_db = 20 * np.log10(noise_estimate)
    threshold_db = noise_estimate_db + x_db

    return 10 ** (threshold_db / 20)


def find_targets(radar_image, min_distance, lobe_width_X, lobe_width_Y, discr_param, noise_threshold=None):
    """Обнаружение целей с улучшенным порогом."""
    if noise_threshold is None:
        noise_threshold = calculate_noise_threshold(radar_image)

    # Создаем маску локальных максимумов
    local_max = ndimage.maximum_filter(radar_image, size=min_distance) == radar_image

    # Применяем порог
    above_threshold = radar_image > noise_threshold

    # Дополнительный фильтр: отбрасываем слишком слабые цели
    strong_targets = radar_image > (noise_threshold*2)

    detected = local_max & above_threshold & strong_targets

    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))

    print(f"Найдено кандидатов: {len(peaks_coords)}")

    return peaks_coords


###############################################################################
# МОДУЛЬ 3: ВЫДЕЛЕНИЕ ОКОН И СЕЧЕНИЙ
###############################################################################

def extract_target_window(radar_image, center_yx, window_size):
    """
    Выделяет окно вокруг цели.
    """
    center_y, center_x = center_yx
    win_h, win_w = window_size
    half_h, half_w = win_h // 2, win_w // 2

    # Определяем границы
    y_start = max(0, center_y - half_h)
    y_end = min(radar_image.shape[0], center_y + half_h)
    x_start = max(0, center_x - half_w)
    x_end = min(radar_image.shape[1], center_x + half_w)

# Вырезаем окно
    window = radar_image[y_start:y_end, x_start:x_end]

    return window


def extract_sections(window):
    """
    Извлекает горизонтальное и вертикальное сечения из окна.
    """
    center_y, center_x = window.shape[0] // 2, window.shape[1] // 2

    horizontal_section = window[center_y, :]  # Центральная строка
    vertical_section = window[:, center_x]  # Центральный столбец

    return horizontal_section, vertical_section


def generate_sinc_signal_from_section(section, window_size, discr_param):
    """
    Преобразует сечение в формат для анализа sinc-функции.
    """
    # Нормируем сечение
    section_linear = np.abs(section)
    section_norm = section_linear / np.max(section_linear)

    # Преобразуем в dB
    section_db = 20 * np.log10(section_norm)

    # Создаем временную ось
    t = np.linspace(-window_size/(2*discr_param), window_size/(2*discr_param), len(section))

    return t, section_db, section_norm

def convert_to_meters(width_samples, direction):
    """Конвертирует ширину из отсчетов в метры."""
    if direction == 'horizontal':
        return width_samples * RANGE_RESOLUTION
    elif direction == 'vertical':
        return width_samples * AZIMUTH_RESOLUTION
    else:
        return width_samples

###############################################################################
# МОДУЛЬ 5: ВИЗУАЛИЗАЦИЯ
###############################################################################

def plot_radar_image(radar_image, detected_peaks, original_targets=None):
    """
    Визуализация радиолокационного изображения с обнаруженными целями.
    """
    plt.figure(figsize=(10, 8))
    plt.imshow(radar_image, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Амплитуда')
    plt.title('Радиолокационное изображение (РЛИ)')
    plt.xlabel('Пиксели (X)')
    plt.ylabel('Пиксели (Y)')

    # Отмечаем обнаруженные цели
    for i, (y, x) in enumerate(detected_peaks):
        plt.plot(x, y, 'bx', markersize=10, markeredgewidth=2,
                 label='Обнаруженные цели' if i == 0 else "")

    # Отмечаем исходные цели (если известны)
    if original_targets:
        for i, (x, y) in enumerate(original_targets):
            plt.plot(x, y, 'go', markersize=8, markeredgewidth=2,
                     label='Исходные цели' if i == 0 else "")

    plt.legend()
    plt.tight_layout()
    plt.show()

#Визуализация анализа цели.
def plot_target_analysis(window, horizontal_section, vertical_section, h_results, v_results, target_id):

    # 1. Окно с целью
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2)
    im = ax1.imshow(window, cmap='hot', interpolation='nearest')
    plt.colorbar(im, ax=ax1)
    ax1.set_title(f'Цель {target_id}: Окно {WINDOW_SIZE}')
    ax1.set_xlabel('X (отн.)')
    ax1.set_ylabel('Y (отн.)')

# Отмечаем сечения
    center_y, center_x = window.shape[0] // 2, window.shape[1] // 2
    ax1.axhline(y=center_y, color='b', linestyle='--', alpha=0.7)
    ax1.axvline(x=center_x, color='g', linestyle='--', alpha=0.7)
    ax1.plot(center_x, center_y, 'wx', markersize=10, markeredgewidth=2)

# 2. Горизонтальное сечение
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=2)
    ax2.plot(horizontal_section, 'g-', linewidth=2, label='Горизонтальное сечение')

    if h_results and h_results['width'] is not None:
        ax2.axvline(x=h_results['wl_norm'], color='r', linestyle=':', alpha=0.7)
        ax2.axvline(x=h_results['wr_norm'], color='r', linestyle=':', alpha=0.7)
        ax2.axhline(y=h_results['threshold_val'], color='r', linestyle='--', alpha=0.7)

    ax2.set_title(f'Горизонтальное сечение')
    ax2.set_xlabel('Отсчёты')
    ax2.set_ylabel('Амплитуда')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

# 3. Вертикальное сечение
    ax3 = plt.subplot2grid((3, 3), (0, 2), rowspan=2)
    ax3.plot(vertical_section, 'b-', linewidth=2, label='Вертикальное сечение')

    if v_results and v_results['width'] is not None:
        ax3.axvline(x=v_results['wl_norm'], color='r', linestyle=':', alpha=0.7)
        ax3.axvline(x=v_results['wr_norm'], color='r', linestyle=':', alpha=0.7)
        ax3.axhline(y=v_results['threshold_val'], color='r', linestyle='--', alpha=0.7)

    ax3.set_title(f'Вертикальное сечение')
    ax3.set_xlabel('Отсчёты')
    ax3.set_ylabel('Амплитуда')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

# 4. Результаты анализа
    ax4 = plt.subplot2grid((3, 3), (2, 2))
    ax4.axis('off')

    if h_results and v_results:
        text_str = f'Цель {target_id} - Результаты:\n\n'
        text_str += f'Горизонтальное сечение:\n'
        text_str += f'  Ширина: {h_results["width"]:.4f}\n'
        text_str += f'  УБЛ: {h_results["classical_pslr"]:.2f} дБ\n'
        text_str += f'  Инт. УБЛ: {h_results["integral_pslr"]:.2f} дБ\n\n'
        text_str += f'Вертикальное сечение:\n'
        text_str += f'  Ширина: {v_results["width"]:.4f}\n'
        text_str += f'  УБЛ: {v_results["classical_pslr"]:.2f} дБ\n'
        text_str += f'  Инт. УБЛ: {v_results["integral_pslr"]:.2f} дБ'

        ax4.text(0.1, 0.9, text_str, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.show()


def plot_target_window(window, target_id):
    """Визуализация окна цели в линейном и логарифмическом масштабах."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Линейный масштаб
    im1 = ax1.imshow(window, cmap='hot', interpolation='nearest')
    plt.colorbar(im1, ax=ax1)
    ax1.set_title(f'Цель {target_id} - Линейный масштаб')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')

    # Логарифмический масштаб (дБ)
    window_db = 20 * np.log10(np.abs(window))
    im2 = ax2.imshow(window_db, cmap='hot', interpolation='nearest')
    plt.colorbar(im2, ax=ax2)
    ax2.set_title(f'Цель {target_id} - Логарифмический масштаб (дБ)')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')

    plt.tight_layout()
    plt.show()


###############################################################################
# ОСНОВНАЯ ФУНКЦИЯ
###############################################################################

def main():
    radar_image = generate_radar_image(TARGETS, IMAGE_SIZE, SINC_WIDTH_X, SINC_WIDTH_Y, DISCR_PARAM, noise_level_db=-40)
    print(f"   Количество исходных целей: {len(TARGETS)}")

    threshold= calculate_noise_threshold(radar_image, window_size=(50, 50), x_db=10)

    detected_peaks = find_targets(radar_image,MIN_DISTANCE, SINC_WIDTH_X, SINC_WIDTH_Y, DISCR_PARAM, threshold)
    print(f"   Обнаружено целей: {len(detected_peaks)}")

    # Визуализация РЛИ
    plot_radar_image(radar_image, detected_peaks, TARGETS)


    print("\nИТОГОВЫЙ ОТЧЕТ")

    for i, target_yx in enumerate(detected_peaks):
        print(f"\nЦель {i + 1} (координаты: {target_yx}):")

        # Выделение окна
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)

        # Визуализация окна в двух масштабах
        plot_target_window(window, i + 1)

# Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window)

# Анализ горизонтального сечения
        print(f"\nГоризонтальное сечение цели {i + 1}")
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(
            horizontal_section, WINDOW_SIZE[0], DISCR_PARAM)
        width_h, classical_pslr_h, integral_pslr_h = analiz_sechenia(t_h, h_signal_db)

# Анализ вертикального сечения
        print(f"\nВертикальное сечение цели {i + 1}")
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(
            vertical_section, WINDOW_SIZE[1], DISCR_PARAM)
        width_v, classical_pslr_v, integral_pslr_v = analiz_sechenia(t_v, v_signal_db)

    #Формирование отчета
        typ_filename = os.path.join(result_folder, "ozenca.typ")

        # Генерируем typst-код
        typ_content = f'''#set page(width: auto, height: auto, margin: 1.5cm)
                              #set text(font: "New Computer Modern", size: 12pt)
                              #show heading: set text(weight: "bold")

                              #align(center)[
          #text(size: 24pt, weight: "bold")[Анализ радиолокационного изображения]
        ]

        #block[
          Размер голограммы: {0} пикселей \
          Количество целей: {len(detected_peaks)}
        ]
        '''

        # Добавляем данные по каждой цели
        for i, target_data in enumerate(detected_peaks):
            target_coords = detected_peaks[i]

            typ_content += f'''
        #align(center)[
          #text(size: 18pt, weight: "bold")[Цель №{i + 1}]
        ]

        Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}

        // Параметры отклика от цели №{i + 1}
        #table(
          columns: 2,
          align: center,
          stroke: (x: 0.5pt, y: 0.5pt),
          inset: 5pt,

          [*Параметр*], [*Значение*],

          [Ширина главного лепестка (гориз.)], [{width_h:.4f}],
          [Максимальный УБЛ (гориз.)], [{classical_pslr_h:.2f} дБ],
          [Интегральный УБЛ (гориз.)], [{integral_pslr_h:.2f} дБ],

          [Ширина главного лепестка (верт.)], [{width_v:.4f}],
          [Максимальный УБЛ (верт.)], [{classical_pslr_v:.2f} дБ],
          [Интегральный УБЛ (верт.)], [{integral_pslr_v:.2f} дБ],
        )
        '''

        with open(typ_filename, 'w', encoding='utf-8') as typ_file:
            typ_file.write(typ_content)
        typ_path = os.path.join(result_folder, "ozenca.typ")
        pdf_path = typ_path.replace('.typ', '.pdf')

        # Вызываем компилятор Typst
        typst_path = r'C:\Users\Gisich_AV\Desktop\typst\typst\typst.exe'
        result = subprocess.run([typst_path, 'compile', typ_path, pdf_path])

main()


import numpy as np
from scipy import stats

def generate_radar_image_with_snr(targets, image_size, lobe_width_X, lobe_width_Y, discr_param, snr_db=20):
    """Генерация РЛИ с гауссовским шумом заданного SNR."""
    clean_image = np.zeros(image_size)
    for target in targets:
        x, y = target
        sinc_target = generate_2d_sinc(x, y, image_size, lobe_width_X, lobe_width_Y, discr_param)
        clean_image += sinc_target
    
    signal_power = np.mean(clean_image**2)
    noise_power = signal_power / (10**(snr_db/10))
    noise_std = np.sqrt(noise_power)
    
    gaussian_noise = np.random.normal(0, noise_std, image_size)
    radar_image = clean_image + gaussian_noise
    
    return radar_image, noise_std

def find_noise_region_simple(radar_image, expected_std, tolerance_db=3, num_samples=50):
    """Поиск шумовой области по соответствию стандартного отклонения."""
    h, w = radar_image.shape
    window_size = (32, 32)
    window_h, window_w = window_size
    
    tolerance_linear = 10**(tolerance_db/20)
    
    for _ in range(num_samples):
        y = np.random.randint(0, h - window_h)
        x = np.random.randint(0, w - window_w)
        
        region = radar_image[y:y+window_h, x:x+window_w]
        region_std = np.std(region)
        
        # Проверяем соответствие стандартного отклонения
        if abs(region_std - expected_std) / expected_std <= (tolerance_linear - 1):
            return region_std
    
    # Если не нашли - используем ожидаемое значение
    return expected_std

def calculate_noise_threshold_simple(radar_image, expected_std, x_db=10):
    """Простой расчет порога через SNR."""
    # Ищем шумовую область
    measured_std = find_noise_region_simple(radar_image, expected_std)
    
    # Порог = уровень_шума + X дБ
    noise_level_db = 20 * np.log10(measured_std + 1e-12)
    threshold_db = noise_level_db + x_db
    threshold_linear = 10**(threshold_db/20)
    
    print(f"Уровень шума: {measured_std:.6f} ({noise_level_db:.2f} дБ)")
    print(f"Порог: {threshold_linear:.6f} ({threshold_db:.2f} дБ)")
    
    return threshold_linear

# В основной функции:
def main():
    SNR_DB = 20
    X_DB = 10
    
    # Генерация с шумом
    radar_image, expected_std = generate_radar_image_with_snr(
        TARGETS, IMAGE_SIZE, SINC_WIDTH_X, SINC_WIDTH_Y, DISCR_PARAM, snr_db=SNR_DB)
    
    # Порог
    noise_threshold = calculate_noise_threshold_simple(radar_image, expected_std, x_db=X_DB)
    
    # Дальше твой существующий код...
