СтерЖенёк(Ферритовый), [05.12.2025 1: 19]
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, stats
import json


def select_polygons_interactive(amplitude_image, max_polygons=10):
    """
    Полуавтоматическое выделение полигонов через интерактивный интерфейс.
    Возвращает список полигонов в виде бинарных масок.
    """
    print("Выделите полигоны кликами. Нажмите 'enter' для завершения выделения.")

    polygons_masks = []
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.imshow(amplitude_image, cmap='gray', aspect='auto')
    ax.set_title(f'Выделите до {max_polygons} полигонов. Клик ЛКМ - точки, ПКМ - завершить полигон')

    for poly_idx in range(max_polygons):
        print(f"\nПолигон {poly_idx + 1}/{max_polygons}")
        print("Кликайте ЛКМ для добавления точек полигона")
        print("Клик ПКМ или нажмите 'enter' для завершения текущего полигона")
        print("Нажмите 'q' для завершения всех выделений")

        points = []

        def onclick(event):
            if event.button == 1:  # ЛКМ
                points.append((int(event.ydata), int(event.xdata)))
                ax.plot(event.xdata, event.ydata, 'ro', markersize=4)
                if len(points) > 1:
                    # Соединяем точки линией
                    last_y, last_x = points[-2]
                    ax.plot([last_x, event.xdata], [last_y, event.ydata], 'r-', linewidth=1)
                fig.canvas.draw()
                print(f"Точка добавлена: ({int(event.ydata)}, {int(event.xdata)})")
            elif event.button == 3:  # ПКМ
                if len(points) >= 3:
                    plt.disconnect(cid)
                    print(f"Полигон {poly_idx + 1} завершен с {len(points)} точками")
                else:
                    print("Нужно минимум 3 точки для полигона")

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.waitforbuttonpress()  # Ждем нажатия enter или q
        plt.disconnect(cid)

        if not points:
            print("Выделение завершено")
            break

        # Создаем маску для полигона
        from skimage.draw import polygon
        poly_mask = np.zeros_like(amplitude_image, dtype=bool)
        rows, cols = polygon([p[0] for p in points], [p[1] for p in points],
                             shape=amplitude_image.shape)
        poly_mask[rows, cols] = True

        # Проверяем, что полигон не слишком маленький
        if np.sum(poly_mask) < 100:  # Минимум 100 пикселей
            print("Полигон слишком маленький, пропускаем...")
            continue

        polygons_masks.append(poly_mask)

        # Отображаем заполненный полигон
        ax.imshow(poly_mask, alpha=0.3, cmap='Reds')
        fig.canvas.draw()

    plt.close(fig)
    return polygons_masks


def compute_polygon_metrics(amplitude_image, complex_image, polygon_mask):
    """
    Вычисление всех возможных метрик для одного полигона.
    """
    # Извлекаем данные только внутри полигона
    mask_indices = np.where(polygon_mask)
    amp_values = amplitude_image[mask_indices]
    complex_values = complex_image[mask_indices]

    if len(amp_values) == 0:
        return None

    # 1. Базовые статистики амплитуды
    amp_mean = np.mean(amp_values)
    amp_std = np.std(amp_values)
    amp_var = np.var(amp_values)
    amp_skew = stats.skew(amp_values.flatten())
    amp_kurtosis = stats.kurtosis(amp_values.flatten())

    # 2. Комплексные статистики
    real_values = np.real(complex_values)
    imag_values = np.imag(complex_values)

    # Проверка круговой симметрии (для АБГШ)
    # Cov(I,Q) ≈ 0 и Var(I) ≈ Var(Q)
    cov_iq = np.cov(real_values, imag_values)[0, 1]
    var_i = np.var(real_values)
    var_q = np.var(imag_values)
    circularity_ratio = min(var_i, var_q) / max(var_i, var_q) if max(var_i, var_q) > 0 else 0


СтерЖенёк(Ферритовый), [05.12.2025 1: 19]
# 3. Распределение амплитуды (тест на Релея)
# Параметр масштаба Релея
sigma_rayleigh = np.sqrt(np.sum(amp_values ** 2) / (2 * len(amp_values)))

# K-S тест против распределения Релея
ks_stat, ks_pvalue = stats.kstest(
    amp_values / amp_mean if amp_mean > 0 else amp_values,  # Нормализованные значения
    'rayleigh',
    args=(sigma_rayleigh,)
)

# 4. Распределение фазы (тест на равномерность)
phases = np.angle(complex_values)
phases_norm = (phases + np.pi) % (2 * np.pi)  # Нормализуем к [0, 2π]

# Тест Рэлея на равномерность фазы
# Для равномерного распределения R ≈ 0
R = np.abs(np.mean(np.exp(1j * phases_norm)))

# 5. Пространственные метрики (текстура)
# Вырезаем ограничивающий прямоугольник полигона для анализа текстуры
rows, cols = np.where(polygon_mask)
if len(rows) > 0 and len(cols) > 0:
    rmin, rmax = rows.min(), rows.max()
    cmin, cmax = cols.min(), cols.max()

    # Берем область с небольшим запасом
    r_start = max(0, rmin - 5)
    r_end = min(amplitude_image.shape[0], rmax + 5)
    c_start = max(0, cmin - 5)
    c_end = min(amplitude_image.shape[1], cmax + 5)

    region = amplitude_image[r_start:r_end, c_start:c_end]
    mask_region = polygon_mask[r_start:r_end, c_start:c_end]

    if region.size > 100:  # Достаточный размер для анализа текстуры
        # Нормализуем область
        if np.max(region) > np.min(region):
            region_norm = (region - np.min(region)) / (np.max(region) - np.min(region))
        else:
            region_norm = region

        # Простые метрики текстуры
        # Градиенты (резкость)
        grad_y, grad_x = np.gradient(region_norm)
        gradient_magnitude = np.sqrt(grad_x2 + grad_y2)
        texture_sharpness = np.mean(gradient_magnitude[mask_region[r_start:r_end, c_start:c_end]])

        # Локальная дисперсия (неоднородность)
        local_var = ndimage.generic_filter(region_norm, np.var, size=3)
        texture_heterogeneity = np.mean(local_var[mask_region[r_start:r_end, c_start:c_end]])
    else:
        texture_sharpness = 0
        texture_heterogeneity = 0
else:
    texture_sharpness = 0
    texture_heterogeneity = 0

# 6. Автокорреляционные метрики
# Одномерная автокорреляция (быстрый метод через FFT)
if len(amp_values) > 100:
    # Берем случайную выборку для скорости
    sample_size = min(1000, len(amp_values))
    sample_idx = np.random.choice(len(amp_values), sample_size, replace=False)
    amp_sample = amp_values[sample_idx]

    # Нормируем
    amp_norm = (amp_sample - np.mean(amp_sample)) / (np.std(amp_sample) + 1e-10)

    # Вычисляем автокорреляцию через FFT
    autocorr = np.correlate(amp_norm, amp_norm, mode='full')
    autocorr = autocorr[len(autocorr) // 2:]  # Берем только положительные лаги
    autocorr = autocorr / autocorr[0]  # Нормируем

    # Характеристики автокорреляции
    # Ширина на уровне 0.5
    try:
        autocorr_width_idx = np.where(autocorr < 0.5)[0][0]
        autocorr_width = autocorr_width_idx
    except:
        autocorr_width = len(autocorr)

    # Энтропия автокорреляции
    autocorr_norm = autocorr / np.sum(autocorr)
    autocorr_entropy = -np.sum(autocorr_norm * np.log(autocorr_norm + 1e-10))
else:
    autocorr_width = 0
    autocorr_entropy = 0

# 7. Классификационные метрики
# Коэффициент вариации
cv = amp_std / amp_mean if amp_mean > 0 else 0

# Отношение среднего к медиане (для выявления выбросов)

СтерЖенёк(Ферритовый), [05.12.2025 1: 19]
median_ratio = amp_mean / np.median(amp_values) if np.median(amp_values) > 0 else 1

metrics = {
    # Базовые статистики
    'pixel_count': len(amp_values),
    'area_pixels': np.sum(polygon_mask),
    'amplitude_mean': float(amp_mean),
    'amplitude_std': float(amp_std),
    'amplitude_var': float(amp_var),
    'amplitude_skew': float(amp_skew),
    'amplitude_kurtosis': float(amp_kurtosis),

    # Комплексные свойства
    'circularity_ratio': float(circularity_ratio),
    'cov_iq': float(cov_iq),
    'var_i': float(var_i),
    'var_q': float(var_q),

    # Тесты распределений
    'rayleigh_sigma': float(sigma_rayleigh),
    'ks_test_pvalue': float(ks_pvalue),
    'phase_uniformity_R': float(R),

    # Пространственные метрики
    'texture_sharpness': float(texture_sharpness),
    'texture_heterogeneity': float(texture_heterogeneity),

    # Автокорреляция
    'autocorr_width': int(autocorr_width),
    'autocorr_entropy': float(autocorr_entropy),

    # Классификационные
    'coefficient_of_variation': float(cv),
    'mean_median_ratio': float(median_ratio),

    # Производные метрики для классификации
    'is_low_variance': float(cv < 0.1),
    'is_high_circularity': float(circularity_ratio > 0.9),
    'is_rayleigh_like': float(ks_pvalue > 0.05),  # Не отвергаем гипотезу о распределении Релея
    'is_uniform_phase': float(R < 0.1),  # R близко к 0 для равномерного распределения
}

return metrics


def analyze_all_polygons(radar_image_complex, polygons_masks):
    """
    Анализ всех полигонов и сохранение результатов.
    """
    amplitude_image = np.abs(radar_image_complex)

    results = {}
    for i, mask in enumerate(polygons_masks):
        print(f"Анализ полигона {i + 1}...", end='\r')
        metrics = compute_polygon_metrics(amplitude_image, radar_image_complex, mask)
        if metrics:
            results[f'polygon_{i + 1}'] = metrics

            # Быстрая классификация по метрикам
            if metrics['is_rayleigh_like'] and metrics['is_uniform_phase'] and metrics['is_high_circularity']:
                results[f'polygon_{i + 1}']['predicted_class'] = 'Шум (АБГШ)'
            elif metrics['coefficient_of_variation'] > 0.3 and metrics['texture_heterogeneity'] > 0.1:
                results[f'polygon_{i + 1}']['predicted_class'] = 'Фон (растительность, город)'
            elif metrics['amplitude_mean'] > 2 * np.mean(amplitude_image):
                results[f'polygon_{i + 1}']['predicted_class'] = 'Цель/Отражение'
            else:
                results[f'polygon_{i + 1}']['predicted_class'] = 'Неизвестно'

    print(f"\nАнализ завершен. Обработано {len(results)} полигонов.")
    return results


# Основная функция для использования
def run_polygon_analysis(hdf5_path, json_path):
    """
    Полный цикл: загрузка, выделение полигонов, анализ.
    """
    # Загрузка данных (используем твои существующие функции)
    radar_image_complex, radar_image = load_radar_image_from_hdf5(hdf5_path)

    # Полуавтоматическое выделение полигонов
    print("Начало выделения полигонов...")
    polygons_masks = select_polygons_interactive(radar_image)

    if not polygons_masks:
        print("Полигоны не выделены. Использую автоматическое разбиение на сетку.")
        # Автоматическое разбиение на 9 равных областей
        h, w = radar_image.shape
        polygons_masks = []
        for i in range(3):
            for j in range(3):
                mask = np.zeros_like(radar_image, dtype=bool)
                y_start = i * h // 3
                y_end = (i + 1) * h // 3
                x_start = j * w // 3
                x_end = (j + 1) * w // 3


СтерЖенёк(Ферритовый), [05.12.2025 1: 19]
mask[y_start:y_end, x_start:x_end] = True
polygons_masks.append(mask)

# Анализ всех полигонов
results = analyze_all_polygons(radar_image_complex, polygons_masks)

# Сохранение результатов
output_file = 'polygon_analysis_results.json'
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False, default=float)

print(f"Результаты сохранены в {output_file}")

# Быстрый отчет
print("\n=== КРАТКИЙ ОТЧЕТ ===")
for poly_name, metrics in results.items():
    print(f"\n{poly_name}: {metrics['predicted_class']}")
    print(f"  Площадь: {metrics['area_pixels']} пикселей")
    print(f"  Средняя амплитуда: {metrics['amplitude_mean']:.2f}")
    print(f"  Коэф. вариации: {metrics['coefficient_of_variation']:.3f}")
    print(f"  R-тест фазы: {metrics['phase_uniformity_R']:.3f} "
          f"{'(равномерно)' if metrics['is_uniform_phase'] else '(не равномерно)'}")
    print(f"  K-S p-value: {metrics['ks_test_pvalue']:.3f} "
          f"{'(Релей)' if metrics['is_rayleigh_like'] else '(не Релей)'}")

return results


# Пример использования в твоем main():
def main():
    # ... существующий код загрузки ...

    # Анализ полигонов
    polygon_results = run_polygon_analysis(HDF5_FILE_PATH, JSON_FILE_PATH)

    # Используем результаты для улучшенного поиска шума
    # Находим полигон, классифицированный как шум
    noise_polygons = []
    for poly_name, metrics in polygon_results.items():
        if metrics.get('predicted_class') == 'Шум (АБГШ)':
            # Получаем координаты полигона
            # (здесь нужна доработка для получения маски из результатов)
            pass

    # ... остальной код ...
############
# РАДИОМЕТРИЧЕСКИЙ АНАЛИЗ - ИСПРАВЛЕННЫЙ ФУНКЦИОНАЛ
############

def calculate_radiometric_sensitivity(targets_energy, noise_power, range_res, azimuth_res, sigma_ref):
    """
    Оценка радиометрической чувствительности для всех целей со статистикой
    """
    nesz_values = []

    for target_energy in targets_energy:
        # Индивидуальный калибровочный коэффициент для каждой цели
        individual_calibration = (target_energy * azimuth_res * range_res) / sigma_ref
        nesz_linear = noise_power / individual_calibration
        nesz_db = 10 * np.log10(nesz_linear)
        nesz_values.append(nesz_db)

    # Статистика
    nesz_mean = np.mean(nesz_values)
    nesz_std = np.std(nesz_values)
    nesz_min = np.min(nesz_values)
    nesz_max = np.max(nesz_values)

    return nesz_values, nesz_mean, nesz_std, nesz_min, nesz_max


def calculate_targets_epr(targets_energy, range_res, azimuth_res, sigma_ref):
    """
    Оценка ЭПР целей через среднюю энергию всех целей
    """
    # Калибровочный коэффициент по средней энергии всех целей
    mean_energy = np.mean(targets_energy)
    calibration_constant = (mean_energy * azimuth_res * range_res) / sigma_ref

    # Расчет ЭПР для каждой цели
    rcs_values = []
    for target_energy in targets_energy:
        rcs_linear = (target_energy * azimuth_res * range_res) / calibration_constant
        rcs_db = 10 * np.log10(rcs_linear) if rcs_linear > 0 else -80
        rcs_values.append(rcs_db)

    return rcs_values, calibration_constant


def calculate_target_energy(target_window_complex, noise_power):
    """
    Расчет энергии цели в главном лепестке с проверкой на фон/шум
    """
    window_real = np.abs(target_window_complex)
    max_amplitude = np.max(window_real)
    threshold = max_amplitude * 10 ** (-3 / 20)
    mask = window_real >= threshold

    # Суммируем энергию в главном лепестке
    target_energy = np.sum(window_real[mask] ** 2)

    # Вычитаем шумовую составляющую
    noise_energy = noise_power * np.sum(mask)
    signal_energy = target_energy - noise_energy

    # Проверка: если сигнал меньше шума - это скорее шум
    if signal_energy <= 0:
        signal_energy = 0

    return signal_energy


############
# ОБНОВЛЕННАЯ ОСНОВНАЯ ФУНКЦИЯ
############

def main():
    # Генерация с шумом
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    detected_peaks = find_targets(radar_image_complex, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)
    Noise_power = calculate_noise_threshold(radar_image_complex)

    # Получаем разрешения (это и есть d_r и d_az)
    range_res, azimuth_res = resolution(T_synth, F_r_discr, Full_velocity, radar_image)

    # Константы
    SIGMA_REF = 268.5  # ЭПР эталонного отражателя

    # Сбор энергий всех целей
    targets_energy = []
    for target_yx in detected_peaks:
        window = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        target_energy = calculate_target_energy(window, Noise_power)
        targets_energy.append(target_energy)

    # 1. ОЦЕНКА РАДИОМЕТРИЧЕСКОЙ ЧУВСТВИТЕЛЬНОСТИ
    nesz_values, nesz_mean, nesz_std, nesz_min, nesz_max = calculate_radiometric_sensitivity(
        targets_energy, Noise_power, range_res, azimuth_res, SIGMA_REF
    )

    # 2. ОЦЕНКА ЭПР ЦЕЛЕЙ  
    rcs_values, calibration_constant = calculate_targets_epr(
        targets_energy, range_res, azimuth_res, SIGMA_REF
    )

    targets_data = []

    for i, target_yx in enumerate(detected_peaks):
        # Выделение окна
        window = extract_target_window(radar_image, target_yx, WINDOW_SIZE)

        # Расчет уровня сигнал/шум
        snr_db = calculate_snr_for_target(window, Noise_power)

        # Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window)


СтерЖенёк(Ферритовый), [01.12.2025 2: 49]
# Анализ горизонтального сечения
t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)

# Анализ вертикального сечения
t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)

# Формируем данные для отчета
target_data = {
    'window_linear': window,
    'snr_db': snr_db,
    'rcs_db': rcs_values[i],
    'nesz_db': nesz_values[i],
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
    f.write(f"Название голограммы: {HOLOGRAM_NAME}\n")
    f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
    f.write(f"Количество целей: {len(detected_peaks)}\n")
    f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
    f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
    f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")
    f.write(f"Расчетная ЭПР эталона (σ_ref): {SIGMA_REF:.2f} м²\n")
    f.write(f"Калибровочный коэффициент: {calibration_constant:.2f}\n")
    f.write(f"Радиометрическая чувствительность (NESZ): {nesz_mean:.2f} дБ\n")
    f.write(f"Мин/макс чувствительность: {nesz_min:.2f}/{nesz_max:.2f} дБ\n")
    f.write(f"СКО чувствительности: {nesz_std:.2f} дБ\n")

    # 1. Сохраняем основное РЛИ с целями
plt.figure(figsize=(9, 8))
plt.imshow(radar_image, cmap='hot', interpolation='nearest')
plt.colorbar(label='Амплитуда')
plt.title('Радиолокационное изображение с обнаруженными целями')
plt.xlabel('Отсчёты по дальности')
plt.ylabel('Отсчёты по азимуту')

# Отмечаем обнаруженные цели
for i, (y, x) in enumerate(detected_peaks):
    plt.plot(x, y, 's', markersize=14, markeredgewidth=2, markeredgecolor='red', markerfacecolor='none',
             linestyle='none', label='Цели' if i == 0 else "")

plt.legend()
radar_image_path = os.path.join(result_folder, "radar_image.png")
plt.savefig(radar_image_path, dpi=150, bbox_inches="tight")
plt.close()

# ВИЗУАЛИЗАЦИЯ РАДИОМЕТРИЧЕСКОЙ ЧУВСТВИТЕЛЬНОСТИ (аналогично вашим графикам)
if nesz_values:
    plt.figure(figsize=(12, 6))
    target_numbers = range(1, len(nesz_values) + 1)

    plt.plot(target_numbers, nesz_values, 'bo-', linewidth=2, markersize=6,
             label=f'Чувствительность (средняя: {nesz_mean:.2f} дБ)')
    plt.axhline(y=nesz_mean, color='red', linestyle='--', linewidth=2,
                label=f'Среднее: {nesz_mean:.2f} дБ')

    # Область ±1 СКО
    plt.fill_between(target_numbers, nesz_mean - nesz_std, nesz_mean + nesz_std,

                     СтерЖенёк(Ферритовый), [01.12.2025 2: 49]
    alpha = 0.2, color = 'red', label = f'±1 СКО ({nesz_std:.2f} дБ)')

    plt.xlabel('Номер цели', fontsize=12)
    plt.ylabel('Радиометрическая чувствительность (дБ)', fontsize=12)
    plt.title('Радиометрическая чувствительность по целям', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Статистика в углу
    stats_text = f'Статистика:\nСреднее: {nesz_mean:.2f} дБ\nМин: {nesz_min:.2f} дБ\nМакс: {nesz_max:.2f} дБ\nСКО: {nesz_std:.2f} дБ'
    plt.annotate(stats_text, xy=(0.02, 0.98), xycoords='axes fraction',
                 verticalalignment='top', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    nesz_plot_path = os.path.join(result_folder, "nesz_plot.png")
    plt.savefig(nesz_plot_path, dpi=150, bbox_inches='tight')
    plt.close()

# 2. Сохраняем данные и графики для каждой цели
for i, target_data in enumerate(targets_data):
    target_id = i + 1
    target_coords = detected_peaks[i]

    # Создаем папку для цели
    target_folder = os.path.join(result_folder, f"target_{target_id}")
    os.makedirs(target_folder, exist_ok=True)

    # Рассчитываем ширину в метрах
    h_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['h_width'],
                                       'horizontal')
    v_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['v_width'],
                                       'vertical')

    # Сохраняем параметры цели
    with open(os.path.join(target_folder, f"target_{target_id}_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Цель №{target_id}\n")
        f.write(f"Положение: азимут {target_coords[0]}, дальность {target_coords[1]}\n\n")
        f.write(f"SNR: {target_data['snr_db']:.2f} дБ\n")
        f.write(f"ЭПР: {target_data['rcs_db']:.2f} дБ м²\n")
        f.write(f"Чувствительность: {target_data['nesz_db']:.2f} дБ\n")

        f.write("Сечение по дальности:\n")
        f.write(f"Ширина главного лепестка: {target_data['h_width']:.4f} отсч. ({h_width_meters:.4f} м)\n")
        f.write(f"Максимальный УБЛ: {target_data['h_pslr']:.2f} дБ\n")
        f.write(f"Интегральный УБЛ: {target_data['h_i_pslr']:.2f} дБ\n\n")

        f.write("Сечение по азимуту:\n")
        f.write(f"Ширина главного лепестка: {target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)\n")
        f.write(f"Максимальный УБЛ: {target_data['v_pslr']:.2f} дБ\n")
        f.write(f"Интегральный УБЛ: {target_data['v_i_pslr']:.2f} дБ\n")

    # Сохраняем массивы данных
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_window_linear.txt"), target_data['window_linear'],
               fmt='%6.6f')
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_horizontal_section.txt"),
               np.column_stack([target_data['h_t'], target_data['h_signal_db']]), fmt="%.6f")
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_vertical_section.txt"),
               np.column_stack([target_data['v_t'], target_data['v_signal_db']]), fmt="%.6f")

    # Окно в линейном масштабе
    plt.figure(figsize=(6, 4))
    plt.imshow(target_data['window_linear'], cmap='hot', interpolation='nearest')
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

СтерЖенёк(Ферритовый), [01.12.2025 2: 49]
plt.colorbar(label='Амплитуда (дБ)')
plt.title(f'Цель {target_id} - Логарифмический масштаб')
plt.xlabel('Дальность')
plt.ylabel('Азимут')
db_path = os.path.join(target_folder, f"target_{target_id}_db.png")
plt.savefig(db_path, dpi=150, bbox_inches='tight')
plt.close()

# Графики сечений (ваш существующий код)
# ... [ваш код графиков сечений] ...

# Формирование Typst отчета
typ_content = f"""
#set_page(width: auto, height: auto, margin: 1.5cm)
#set text(font: "New Computer Modern", size: 12pt, lang: "ru")
#show heading: set text(weight: "bold")

#align(center)[
#text(size: 24pt, weight: "bold")[Анализ радиолокационного изображения]
]

#align(center)[
#text(size: 12pt)[Голограмма: {HOLOGRAM_NAME}]
]

#align(center)[
#text(size: 12pt)[Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей, количество целей: {len(detected_peaks)}]
]

#align(center)[
#text(size: 12pt)[Радиометрическая чувствительность: средняя {nesz_mean:.2f} дБ, СКО {nesz_std:.2f} дБ]
]

#align(center)[
#figure(
    image("radar_image.png"),
    caption: [Радиолокационное изображение с обнаруженными целями]
)
"""

if nesz_values:
    typ_content += f"""
#align(center)[
#figure(
    image("nesz_plot.png"),
    caption: [Радиометрическая чувствительность по целям. Среднее: {nesz_mean:.2f} дБ, СКО: {nesz_std:.2f} дБ]
)
]
"""

# Добавляем данные по каждой цели
for i, target_data in enumerate(targets_data):
    target_id = i + 1
    target_coords = detected_peaks[i]

    h_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['h_width'],
                                       'horizontal')
    v_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['v_width'],
                                       'vertical')

    typ_content += f"""
#pagebreak()
#align(center)[
#text(size: 18pt, weight: "bold")[Цель №{target_id}]
]

#align(center)[
#text(size: 12pt)[Положение: азимут: {target_coords[0]}, дальность: {target_coords[1]}]
]

#align(center)[
#text(size: 12pt)[SNR: {target_data['snr_db']:.2f} дБ, ЭПР: {target_data['rcs_db']:.2f} дБ м², Чувствительность: {target_data['nesz_db']:.2f} дБ]
]

// ... [остальная часть вашего Typst кода для целей] ...
"""

# Сохранение и компиляция Typst
typ_filename = os.path.join(result_folder, "ozenca.typ")
with open(typ_filename, 'w', encoding='utf-8') as typ_file:
    typ_file.write(typ_content)

pdf_path = typ_filename.replace('.typ', '.pdf')

# Компиляция в PDF
try:
    typst_path = r'C:\Users\Babaylova\Documents\typst\typst.exe'
    subprocess.run([typst_path, 'compile', typ_filename, pdf_path], capture_output=True, check=True)
except Exception as e:
    pass

if name == "main":
    main()

СтерЖенёк(Ферритовый), [01.12.2025 2: 53]

############
# РАДИОМЕТРИЧЕСКАЯ КОРРЕКЦИЯ ПО МЕТОДУ ICEYE
############

def apply_radiometric_correction(radar_image_complex, json_data, platform_type='airborne'):
    """
    Применяет радиометрическую коррекцию с поправкой на угол падения
    по аналогии с алгоритмом ICEYE

    Parameters:
    - radar_image_complex: комплексное РЛИ
    - json_data: данные из JSON файла
    - platform_type: 'airborne' (авиа) или 'spaceborne' (космос)
    """

    # Получаем геометрические параметры из JSON
    try:
        # Время задержки для первого и последнего отсчета по дальности
        first_delay = json_data['image_description']['first_range_sample_delay']
        last_delay = json_data['image_description']['last_range_sample_delay']

        # Высота платформы (для авиационных данных берем из ellipsoid_params)
        if platform_type == 'airborne':
            platform_height = json_data['ellipsoid_params']['local_height']
        else:
            # Для космических данных может быть другой параметр
            platform_height = json_data.get('platform_height', 500000)  # по умолчанию 500 км

        # Геометрические параметры эллипсоида
        major_axis = json_data['ellipsoid_params']['major_semiaxis']
        minor_axis = json_data['ellipsoid_params']['minor_semiaxis']

    except KeyError as e:
        print(f"Отсутствует необходимый параметр в JSON: {e}")
        return radar_image_complex

    # Размеры изображения
    height, width = radar_image_complex.shape

    # Создаем матрицу углов падения
    incidence_angles = np.zeros((height, width))

    # Расчет углов падения для каждого пикселя
    for i in range(height):
        for j in range(width):
            # Расчет наклонной дальности для текущего пикселя
            range_delay = first_delay + (last_delay - first_delay) * (j / width)
            slant_range = range_delay * 3e8 / 2  # преобразуем время в расстояние

            # Расчет угла падения (упрощенная геометрическая модель)
            # Для авиационных данных
            if platform_type == 'airborne':
                # Упрощенный расчет для плоской Земли
                ground_range = np.sqrt(slant_range2 - platform_height2)
                incidence_angle = np.arcsin(platform_height / slant_range)
            else:
                # Для космических данных с учетом кривизны Земли
                Re = major_axis  # радиус Земли
                sat_radius = Re + platform_height
                # Расчет по теореме косинусов
                cos_theta = (sat_radius2 + slant_range2 - Re ** 2) / (2 * sat_radius * slant_range)
                incidence_angle = np.arccos(cos_theta)

            incidence_angles[i, j] = np.degrees(incidence_angle)  # переводим в градусы

    # Применяем радиометрическую коррекцию по формуле ICEYE
    # Мощность исходного сигнала
    power_original = np.abs(radar_image_complex) ** 2

    # Поправочный коэффициент на основе угла падения
    correction_factor = np.sin(np.radians(incidence_angles))

    # Избегаем деления на ноль
    correction_factor = np.where(correction_factor < 0.01, 0.01, correction_factor)

    # Применяем коррекцию к мощности
    power_corrected = power_original / correction_factor

    # Восстанавливаем комплексное изображение с сохранением фазы
    phase_original = np.angle(radar_image_complex)
    amplitude_corrected = np.sqrt(power_corrected)
    complex_corrected = amplitude_corrected * np.exp(1j * phase_original)

    # Визуализация углов падения для отладки
    plt.figure(figsize=(10, 6))
    plt.imshow(incidence_angles, cmap='jet', aspect='auto')
    plt.colorbar(label='Угол падения (градусы)')
    plt.title('Матрица углов падения')
    plt.savefig(os.path.join(result_folder, "incidence_angles.png"), dpi=150, bbox_inches='tight')


СтерЖенёк(Ферритовый), [01.12.2025 2: 53]
plt.close()

# Визуализация поправочного коэффициента
plt.figure(figsize=(10, 6))
plt.imshow(correction_factor, cmap='jet', aspect='auto')
plt.colorbar(label='Поправочный коэффициент')
plt.title('Матрица поправочных коэффициентов')
plt.savefig(os.path.join(result_folder, "correction_factors.png"), dpi=150, bbox_inches='tight')
plt.close()

return complex_corrected, incidence_angles, correction_factor


############
# ОБНОВЛЕННАЯ ОСНОВНАЯ ФУНКЦИЯ
############

def main():
    # Загрузка данных
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)

    # Загрузка JSON параметров
    with open(JSON_FILE_PATH, 'r', encoding='utf-8') as f:
        json_data = json.load(f)

    # ПРИМЕНЕНИЕ РАДИОМЕТРИЧЕСКОЙ КОРРЕКЦИИ
    radar_image_corrected, incidence_angles, correction_factors = apply_radiometric_correction(
        radar_image_complex, json_data, platform_type='airborne'
    )

    # Обновляем изображение для дальнейшего анализа
    radar_image_complex = radar_image_corrected
    radar_image = np.abs(radar_image_complex)

    # Продолжаем основной анализ...
    detected_peaks = find_targets(radar_image_complex, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)
    Noise_power = calculate_noise_threshold(radar_image_complex)

    # ... остальная часть вашего кода без изменений ...

    # ДОБАВЛЯЕМ ИНФОРМАЦИЮ О КОРРЕКЦИИ В ОТЧЕТ
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Название голограммы: {HOLOGRAM_NAME}\n")
        f.write(f"Применена радиометрическая коррекция: ДА\n")
        f.write(f"Метод коррекции: ICEYE-style с поправкой на угол падения\n")
        f.write(f"Диапазон углов падения: {np.min(incidence_angles):.1f} - {np.max(incidence_angles):.1f} градусов\n")
        f.write(
            f"Диапазон поправочных коэффициентов: {np.min(correction_factors):.3f} - {np.max(correction_factors):.3f}\n")
        # ... остальные параметры ...

    # ДОБАВЛЯЕМ ГРАФИКИ КОРРЕКЦИИ В ОТЧЕТ
    typ_content += f"""
#align(center)[
#text(size: 16pt, weight: "bold")[Результаты радиометрической коррекции]
]

#align(center)[
#table(
    columns: (1.5cm, 1.5cm),
    align: center,
    stroke: (x: 0.2pt, y: 0.2pt),
    inset: 5pt,
    [
        #figure(
            image("incidence_angles.png", width: 100%),
            caption: [Матрица углов падения]
        )
    ],
    [
        #figure(
            image("correction_factors.png", width: 100%),
            caption: [Поправочные коэффициенты]
        )
    ]
]
"""


# УПРОЩЕННАЯ ВЕРСИЯ ДЛЯ БЫСТРОГО ТЕСТИРОВАНИЯ
def quick_radiometric_correction(radar_image_complex, incidence_angle_deg=45):
    """
    Упрощенная версия коррекции с постоянным углом падения
    """
    # Преобразуем угол в радианы
    theta_rad = np.radians(incidence_angle_deg)

    # Применяем коррекцию по формуле ICEYE
    power_original = np.abs(radar_image_complex) ** 2
    power_corrected = power_original / np.sin(theta_rad)

    # Сохраняем фазу
    phase_original = np.angle(radar_image_complex)
    amplitude_corrected = np.sqrt(power_corrected)
    complex_corrected = amplitude_corrected * np.exp(1j * phase_original)

    return complex_corrected


СтерЖенёк(Ферритовый), [01.12.2025 3: 19]

def apply_radiometric_correction(radar_image_complex, json_data):
    # 1. Извлекаем параметры из JSON
    first_delay = json_data['image_description']['first_range_sample_delay']
    last_delay = json_data['image_description']['last_range_sample_delay']
    platform_height = json_data['ellipsoid_params']['local_height']

    # 2. Подготавливаем матрицы
    height, width = radar_image_complex.shape
    incidence_angles = np.zeros((height, width))

    # 3. Расчет для каждого пикселя
    for i in range(height):
        for j in range(width):
            # Расчет наклонной дальности
            range_delay = first_delay + (last_delay - first_delay) * (j / width)
            slant_range = range_delay * 3e8 / 2

            # Расчет угла падения (авиационная модель)
            incidence_angle = np.arcsin(platform_height / slant_range)
            incidence_angles[i, j] = np.degrees(incidence_angle)

    # 4. Применение коррекции
    power_original = np.abs(radar_image_complex) ** 2
    correction_factor = np.sin(np.radians(incidence_angles))
    power_corrected = power_original / correction_factor

    # 5. Восстановление комплексного изображения
    phase_original = np.angle(radar_image_complex)
    amplitude_corrected = np.sqrt(power_corrected)
    complex_corrected = amplitude_corrected * np.exp(1j * phase_original)

    return complex_corrected, incidence_angles, correction_factor


СтерЖенёк(Ферритовый), [01.12.2025 3: 31]

############
# УПРОЩЕННАЯ ПРОВЕРКА ШУМ/ФОН
############

def check_noise_vs_background(region_complex, noise_threshold=0.1):
    """
    Простая проверка: шум или фон на основе математического ожидания

    Parameters:
    - region_complex: проверяемая область (комплексная)
    - noise_threshold: порог (доля от средней амплитуды)

    Returns:
    - bool: True если шум, False если фон
    """

    # Вычисляем математическое ожидание реальной и мнимой частей
    mean_real = np.mean(np.real(region_complex))
    mean_imag = np.mean(np.imag(region_complex))
    mean_amplitude = np.mean(np.abs(region_complex))

    # Простая проверка: если средние близки к нулю - это шум
    # Порог: среднее меньше чем X% от средней амплитуды
    threshold_value = noise_threshold * mean_amplitude

    is_noise = (abs(mean_real) < threshold_value and
                abs(mean_imag) < threshold_value)

    return is_noise


def calculate_noise_threshold(radar_image_complex):
    """Расчет порога на основе мощности шума - УПРОЩЕННАЯ ВЕРСИЯ"""

    # Находим шумовую область (ваша существующая функция)
    noise_window, amplitudes, phases, match_quality = find_rayleigh_uniform_region(radar_image_complex)

    # Простая проверка что это шум
    is_confirmed_noise = check_noise_vs_background(noise_window)

    # Если это не шум, берем минимальную по амплитуде область
    if not is_confirmed_noise:
        h, w = radar_image_complex.shape
        min_amplitude = float('inf')
        best_region = noise_window

        for attempt in range(5):  # всего 5 попыток
            y = np.random.randint(0, h - 50)
            x = np.random.randint(0, w - 50)
            test_region = radar_image_complex[y:y + 50, x:x + 50]
            test_amplitude = np.mean(np.abs(test_region))

            if test_amplitude < min_amplitude:
                min_amplitude = test_amplitude
                best_region = test_region

        noise_window = best_region
        amplitudes = np.abs(noise_window).flatten()

    # Вычисляем мощность шума
    noise_power = np.mean(amplitudes ** 2)

    return noise_power


############
# ОБНОВЛЕННАЯ ОСНОВНАЯ ФУНКЦИЯ (минимальные изменения)
############

def main():
    # Загрузка данных (ваш существующий код)
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)

    # Расчет порога шума (просто вызываем функцию)
    Noise_power = calculate_noise_threshold(radar_image_complex)

    # Обнаружение целей (ваш существующий код)
    detected_peaks = find_targets(radar_image_complex, MIN_DISTANCE, THRESHOLD_OFFSET_DB)

    # Анализ каждой цели
    targets_data = []
    for i, target_yx in enumerate(detected_peaks):
        window = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)

        # ПРОСТАЯ ПРОВЕРКА: шум или фон вокруг цели
        is_background = not check_noise_vs_background(window)

        # Расчет SNR (ваш существующий код)
        snr_db = calculate_snr_for_target(window, Noise_power)

        # Формируем данные цели
        target_data = {
            'window_linear': np.abs(window),
            'snr_db': snr_db,
            'is_on_background': is_background,
            # ... остальные ваши параметры без изменений ...
        }
        targets_data.append(target_data)

    # Просто добавляем одну строку в отчет
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        # ... ваши существующие записи ...
        f.write(f"Мощность шума: {Noise_power:.6f}\n")

        # Простая статистика по целям на фоне
        targets_on_background = sum(1 for t in targets_data if t['is_on_background'])
        f.write(f"Целей на фоне: {targets_on_background} из {len(targets_data)}\n")