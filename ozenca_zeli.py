import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from sinc_interpolation import sinc_interpolation
from find_main_lobe_width import find_main_lobe_width
from calculate_sidelobe_levels import calculate_sidelobe_levels
import matplotlib

# Настройки для качественной визуализации
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['figure.figsize'] = (12, 8)

###############################################################################
# КОНФИГУРАЦИЯ СИСТЕМЫ
###############################################################################

# Априорные параметры РЛИ
IMAGE_SIZE = (512, 512)  # Размер радиолокационного изображения
TARGETS = [(150, 200), (350, 300), (400, 100), (400, 10)]  # Координаты целей (x, y)
SINC_WIDTH_X = 3.0  # Ширина sinc-функции в горизонтальном сечении
SINC_WIDTH_Y = 3.0 # Ширина sinc-функции в вертикальном сечении

# Параметры обнаружения
PEAK_DETECTION = {
    'min_distance': 50,
    'threshold_abs': 0.1,
    'sigma': 1.0
}

# Параметры окон
WINDOW_SIZE = (32, 32)  # Размер окна анализа вокруг цели

F_DISCR = 500
SINC_WIDTH = 3.0
TIME_RANGE = 10.0
KERNEL_SIZE =8
SUBSAMPLES = 16
INTERP_FACTOR = 10
THRESHOLD_DB = -3

###############################################################################
# МОДУЛЬ 1: ГЕНЕРАЦИЯ РАДИОЛОКАЦИОННОГО ИЗОБРАЖЕНИЯ
###############################################################################

def generate_2d_sinc(x0, y0, size, lobe_width_X, lobe_width_Y):
    """
    Генерирует двумерную sinc-функцию как произведение двух одномерных sinc.
    """
    # Создаем координатные сетки
    y, x = np.ogrid[0:size[0], 0:size[1]]

    # Вычисляем нормализованные расстояния от центра
    x_norm = (x - x0) / lobe_width_X
    y_norm = (y - y0) / lobe_width_Y

    # Создаем двумерную sinc как произведение двух одномерных
    sinc_x = np.sinc(x_norm)
    sinc_y = np.sinc(y_norm)
    sinc_2d = sinc_x * sinc_y

    return sinc_2d


def generate_radar_image(targets, image_size, lobe_width_X, lobe_width_Y ):
    """
    Генерирует радиолокационное изображение с целями в виде 2D sinc-функций.
    """
    radar_image = np.zeros(image_size)

    for target in targets:
        x, y = target
        sinc_target = generate_2d_sinc(x, y, image_size, lobe_width_X, lobe_width_Y)
        radar_image += sinc_target

    return radar_image


###############################################################################
# МОДУЛЬ 2: ОБНАРУЖЕНИЕ ЦЕЛЕЙ
###############################################################################

def find_targets(radar_image, min_distance=50, threshold_abs=0.1, sigma=1.0):
    """
    Обнаруживает цели в радиолокационном изображении.
    """
    # Сглаживаем изображение для подавления шумов
    smoothed = ndimage.gaussian_filter(radar_image, sigma=sigma)

    # Создаем маску локальных максимумов
    local_max = ndimage.maximum_filter(smoothed, size=min_distance) == smoothed

    # Применяем порог
    above_threshold = smoothed > threshold_abs
    detected = local_max & above_threshold

    # Получаем координаты обнаруженных целей
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))

    return peaks_coords


###############################################################################
# МОДУЛЬ 3: ВЫДЕЛЕНИЕ ОКОН И СЕЧЕНИЙ
###############################################################################

def extract_target_window(radar_image, center_yx, window_size=(32, 32)):
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

# Дополняем нулями при необходимости
    if window.shape != window_size:
        padded_window = np.zeros(window_size)
        padded_window[:window.shape[0], :window.shape[1]] = window
        return padded_window

    return window


def extract_sections(window):
    """
    Извлекает горизонтальное и вертикальное сечения из окна.
    """
    center_y, center_x = window.shape[0] // 2, window.shape[1] // 2

    horizontal_section = window[center_y, :]  # Центральная строка
    vertical_section = window[:, center_x]  # Центральный столбец

    return horizontal_section, vertical_section

    y_end = min(radar_image.shape[0], center_y + half_h)
    x_start = max(0, center_x - half_w)
    x_end = min(radar_image.shape[1], center_x + half_w)

# Вырезаем окно
    window = radar_image[y_start:y_end, x_start:x_end]

# Дополняем нулями при необходимости
    if window.shape != window_size:
        padded_window = np.zeros(window_size)
        padded_window[:window.shape[0], :window.shape[1]] = window
        return padded_window

    return window


def generate_sinc_signal_from_section(section, time_range=10.0):
    """
    Преобразует сечение в формат для анализа sinc-функции.
    """
    # Нормируем сечение
    section_linear = np.abs(section)
    section_norm = section_linear / np.max(section_linear)

    # Преобразуем в dB
    section_db = 20 * np.log10(section_norm+ 1e-12)

    # Создаем временную ось
    t = np.linspace(-time_range, time_range, len(section))

    return t, section_db, section_norm


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


###############################################################################
# ОСНОВНАЯ ФУНКЦИЯ
###############################################################################

def main():
    radar_image = generate_radar_image(TARGETS, IMAGE_SIZE, SINC_WIDTH_X, SINC_WIDTH_Y)
    print(f"   Количество исходных целей: {len(TARGETS)}")

    detected_peaks = find_targets(radar_image, **PEAK_DETECTION)
    print(f"   Обнаружено целей: {len(detected_peaks)}")

    # 3. Визуализация РЛИ
    plot_radar_image(radar_image, detected_peaks, TARGETS)

    all_results = []

    for i , target_yx in enumerate(detected_peaks):

        # Выделение окна
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)

# Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window)

# Анализ горизонтального сечения
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(
            horizontal_section, TIME_RANGE)
        t_h_interp, h_signal_interp = sinc_interpolation(t_h, h_signal_db, KERNEL_SIZE, SUBSAMPLES, TIME_RANGE, F_DISCR, INTERP_FACTOR)

        wl_h, wr_h, width_h, left_points_h, right_points_h = find_main_lobe_width(
            t_h_interp, h_signal_interp, THRESHOLD_DB)

        classical_pslr_h, integral_pslr_h = calculate_sidelobe_levels(
            t_h_interp, h_signal_interp, SINC_WIDTH_X)

# Нормализуем координаты для отображения в окне
        wl_norm_h = (wl_h + TIME_RANGE) / (2 * TIME_RANGE) * len(horizontal_section)
        wr_norm_h = (wr_h + TIME_RANGE) / (2 * TIME_RANGE) * len(horizontal_section)

        h_results = {
            'width': width_h, 'classical_pslr': classical_pslr_h,
            'integral_pslr': integral_pslr_h, 'wl_norm': wl_norm_h,
            'wr_norm': wr_norm_h, 'threshold_val': 10 ** (THRESHOLD_DB/ 20)
        }

# Анализ вертикального сечения
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(
            vertical_section, TIME_RANGE)
        t_v_interp, v_signal_interp = sinc_interpolation(t_v, v_signal_db, KERNEL_SIZE, SUBSAMPLES, TIME_RANGE, F_DISCR, INTERP_FACTOR)

        wl_v, wr_v, width_v, left_points_v, right_points_v = find_main_lobe_width(
            t_v_interp, v_signal_interp, THRESHOLD_DB)

        classical_pslr_v, integral_pslr_v = calculate_sidelobe_levels(
            t_v_interp, v_signal_interp, SINC_WIDTH_Y)

        wl_norm_v = (wl_v + TIME_RANGE) / (2 * TIME_RANGE) * len(vertical_section)
        wr_norm_v = (wr_v + TIME_RANGE) / (2 * TIME_RANGE) * len(vertical_section)

        v_results = {
            'width': width_v, 'classical_pslr': classical_pslr_v,
            'integral_pslr': integral_pslr_v, 'wl_norm': wl_norm_v,
            'wr_norm': wr_norm_v, 'threshold_val': 10 ** (THRESHOLD_DB / 20)
        }

# Сохраняем результаты
        all_results.append({
        'target_id': i + 1,
        'coordinates': target_yx,
        'horizontal_results': h_results,
        'vertical_results': v_results
        })

# Визуализация анализа цели
        plot_target_analysis(window, horizontal_section, vertical_section, h_results, v_results, i + 1)

    # 5. Итоговый отчет
        print("ИТОГОВЫЙ ОТЧЕТ")

        for result in all_results:
            print(f"\nЦель {result['target_id']} (координаты: {result['coordinates']}):")
            h = result['horizontal_results']
            v = result['vertical_results']

            print(f"  Горизонтальное сечение:")
            print(f"    Ширина лепестка: {h['width']:.4f}")
            print(f"    Классический УБЛ: {h['classical_pslr']:.2f} дБ")
            print(f"    Интегральный УБЛ: {h['integral_pslr']:.2f} дБ")

            print(f"  Вертикальное сечение:")
            print(f"    Ширина лепестка: {v['width']:.4f}")
            print(f"    Классический УБЛ: {v['classical_pslr']:.2f} дБ")
            print(f"    Интегральный УБЛ: {v['integral_pslr']:.2f} дБ")

    return all_results

main()