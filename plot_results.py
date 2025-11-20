import numpy as np

def calculate_noise_threshold(radar_image_complex, x_db):
    """Расчет порога на основе мощности шума с возвратом параметров шума для SNR"""

    # Находим шумовую область
    noise_window, amplitudes, phases, match_quality = find_rayleigh_uniform_region(radar_image_complex)

    # Вычисляем мощность шума
    noise_power = np.mean(amplitudes ** 2)

    # Вычисляем среднеквадратичное значение шума
    noise_rms = np.sqrt(noise_power)

    # Устанавливаем порог
    threshold_db = 20 * np.log10(noise_rms) + x_db
    threshold_linear = 10 ** (threshold_db / 20)

    # === ДОБАВЛЕНО: Возвращаем параметры шума для расчета SNR ===
    return threshold_linear, noise_rms, noise_power, amplitudes

# === ДОБАВЛЕНИЕ РАСЧЕТА SNR ===

def calculate_snr_for_target(target_window_complex, noise_amplitudes):
    """
    Расчет SNR для цели

    Parameters:
    -----------
    target_window_complex : ndarray
        Комплексное окно вокруг цели
    noise_amplitudes : ndarray  
        Амплитуды шумовой области

    Returns:
    --------
    snr_db : float
        SNR в дБ
    """

    # Амплитуды цели
    target_amplitudes = np.abs(target_window_complex).flatten()

    # Мощность сигнала (средняя мощность в окне цели)
    signal_power = np.mean(target_amplitudes ** 2)

    # Мощность шума
    noise_power = np.mean(noise_amplitudes ** 2)

    # SNR в линейной области
    snr_linear = signal_power / noise_power

    # SNR в дБ
    snr_db = 10 * np.log10(snr_linear)

    return snr_db


def calculate_snr_peak(target_window_complex, noise_rms):
    """
    Расчет пикового SNR (отношение максимальной амплитуды цели к RMS шума)
    """
    target_amplitudes = np.abs(target_window_complex)
    peak_amplitude = np.max(target_amplitudes)
    snr_peak_db = 20 * np.log10(peak_amplitude / noise_rms)
    return snr_peak_db


# === ИЗМЕНЕНИЯ В ОСНОВНОЙ ФУНКЦИИ main() ===

def main():
    # Генерация с шумом
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)

    # === ИСПРАВЛЕНИЕ: Используем отладочную версию для проверки фаз ===
    print("=== ПРОВЕРКА КОМПЛЕКСНОГО ИЗОБРАЖЕНИЯ И ФАЗ ===")
    noise_window, noise_amplitudes, noise_phases, match_quality = find_rayleigh_uniform_region_debug(
        radar_image_complex)

    # Расчет порога и параметров шума
    threshold_linear, noise_rms, noise_power = calculate_noise_threshold(radar_image_complex, THRESHOLD_OFFSET_DB)

    detected_peaks = find_targets(radar_image, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)

    targets_data = []

    for i, target_yx in enumerate(detected_peaks):
        # === ИЗМЕНЕНИЕ: Выделяем КОМПЛЕКСНОЕ окно для расчета SNR и фаз ===
        window_complex = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        window_amplitude = np.abs(window_complex)  # Амплитудное окно для визуализации

        # === ДОБАВЛЕНИЕ: Расчет SNR для цели ===
        snr_db = calculate_snr_for_target(window_complex, noise_amplitudes)
        snr_peak_db = calculate_snr_peak(window_complex, noise_rms)

        print(f"Цель {i + 1}: SNR = {snr_db:.2f} дБ, Пиковый SNR = {snr_peak_db:.2f} дБ")

        # Извлечение сечений из амплитудного окна
        horizontal_section, vertical_section = extract_sections(window_amplitude)

        # Анализ горизонтального сечения
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)

        # Анализ вертикального сечения
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)

        # Формируем данные для отчета
        target_data = {
            'window_linear': window_amplitude,  # Используем амплитудное окно для визуализации
            'window_complex': window_complex,  # Сохраняем комплексное для анализа
            'snr_db': snr_db,  # === ДОБАВЛЕНО: SNR ===
            'snr_peak_db': snr_peak_db,  # === ДОБАВЛЕНО: Пиковый SNR ===
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

    # ... остальная часть функции main без изменений до генерации отчета ...

    # === ИЗМЕНЕНИЕ В ОТЧЕТЕ: Добавляем SNR ===
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

        # Сохраняем параметры цели С ДОБАВЛЕННЫМ SNR
        with open(os.path.join(target_folder, f"target_{target_id}_params.txt"), 'w', encoding='utf-8') as f:
            f.write(f"Цель №{target_id}\n")
            f.write(f"Положение: азимут {target_coords[0]}, дальность {target_coords[1]}\n\n")

            f.write("=== ПАРАМЕТРЫ КАЧЕСТВА ===\n")
            f.write(f"SNR: {target_data['snr_db']:.2f} дБ\n")
            f.write(f"Пиковый SNR: {target_data['snr_peak_db']:.2f} дБ\n\n")

            f.write("Сечение по дальности:\n")
            f.write(f"Ширина главного лепестка: {target_data['h_width']:.4f} отсч. ({h_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['h_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['h_i_pslr']:.2f} дБ\n\n")

            f.write("Сечение по азимуту:\n")
            f.write(f"Ширина главного лепестка: {target_data['v_width']:.4f} отсч. ({v_width_meters:.4f} м)\n")
            f.write(f"Максимальный УБЛ: {target_data['v_pslr']:.2f} дБ\n")
            f.write(f"Интегральный УБЛ: {target_data['v_i_pslr']:.2f} дБ\n")

        # ... остальная часть сохранения данных без изменений ...

    # === ИЗМЕНЕНИЕ В TYPST ОТЧЕТЕ: Добавляем SNR ===
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
#text(size: 12pt)[Уровень шума: {20 * np.log10(noise_rms):.2f} дБ, RMS шума: {noise_rms:.6f}]
]

#align(center)[
#figure(
    image("radar_image.png"),
    caption: [Радиолокационное изображение с обнаруженными целями]
)
]
"""

    # Добавляем данные по каждой цели с SNR
    for i, target_data in enumerate(targets_data):
        target_id = i + 1
        target_coords = detected_peaks[i]
        target_folder = f"target_{target_id}"

        # Рассчитываем ширину в метрах для таблицы
        h_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['h_width'],
                                           'horizontal')
        v_width_meters = convert_to_meters(T_synth, F_r_discr, Full_velocity, radar_image, target_data['v_width'],
                                           'vertical')

        typ_content += f"""
#align(center)[
#text(size: 18pt, weight: "bold")[Цель №{target_id}]
]

#align(center)[
#text(size: 12pt)[Положение: азимут {target_coords[0]}, дальность {target_coords[1]}]
]

СтерЖенёк (Ферритовый), [21.11.2025 0:44]
// Параметры качества цели
#align(center)[
#table(
    columns: 2,
    align: center,
    stroke: (x: 0.5pt, y: 0.5pt),
    inset: 5pt,
    [*Параметр*], [*Значение*],
    [SNR], [{target_data['snr_db']:.2f} дБ],
    [Пиковый SNR], [{target_data['snr_peak_db']:.2f} дБ],
)
]

// Визуализация окна цели
#align(center)[
#table(
    columns: (14cm, 14cm), 
    align: center, 
    stroke: (x: 0.2pt, y: 0.2pt), 
    inset: 5pt,
    [
        #figure(
            image("target_{target_id}/target_{target_id}_linear.png", width: 100%),
            caption: [Окно цели - линейный масштаб]
        )
    ],
    [
        #figure(
            image("target_{target_id}/target_{target_id}_db.png", width: 95%),
            caption: [Окно цели - логарифмический масштаб]
        )
    ]
)]

// Сечения цели
#align(center)[
#table(
    columns: (14cm, 14cm), 
    align: center, 
    stroke: (x: 0.2pt, y: 0.2pt),
    inset: 5pt,
    [
        #figure(
            image("target_{target_id}/target_{target_id}_horizontal.png", width: 100%),
            caption: [Сечение по дальности]
        )
    ],
    [
        #figure(
            image("target_{target_id}/target_{target_id}_vertical.png", width: 91%),
            caption: [Сечение по азимуту]
        )
    ]
)]

#align(center)[
#table(
    columns: 2,
    align: center,
    stroke: (x: 0.5pt, y: 0.5pt),
    inset: 5pt,
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
)]
"""

    # ... остальная часть генерации отчета без изменений ...