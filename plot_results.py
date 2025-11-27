
############
# РАДИОМЕТРИЧЕСКАЯ КАЛИБРОВКА И ЧУВСТВИТЕЛЬНОСТЬ
############

def calculate_calibration_coefficient(radar_image_complex, calibration_targets_coords,
                                      sigma_ref, d_az, d_r, window_size=(128, 128)):
    """
    Расчет калибровочного коэффициента по эталонным целям
    """
    calibration_energies = []

    for target_coords in calibration_targets_coords:
        window = extract_target_window(radar_image_complex, target_coords, window_size)
        window_energy = np.sum(np.abs(window) ** 2)
        calibration_energies.append(window_energy)

    # Усредненная энергия эталонных целей
    E_mean_linear = np.mean(calibration_energies)
    E_mean_db = 10 * np.log10(E_mean_linear)

    # Расчет калибровочного коэффициента
    K = E_mean_linear * d_az * d_r / sigma_ref

    return K, E_mean_db


def calculate_nesz(noise_power, K, d_az, d_r):
    """
    Расчет радиометрической чувствительности (NESZ)
    """
    # NESZ в линейных единицах
    nesz_linear = noise_power * d_az * d_r / K
    # Преобразование в дБ
    nesz_db = 10 * np.log10(nesz_linear)

    return nesz_db, nesz_linear


def calculate_sigma0_calibrated(pixel_power, noise_power, K, d_az, d_r):
    """
    Расчет калиброванной УЭПР для пикселя
    """
    # Вычитание шума и калибровка
    sigma0_linear = (pixel_power - noise_power) * d_az * d_r / K
    sigma0_db = 10 * np.log10(sigma0_linear) if sigma0_linear > 0 else -80

    return sigma0_db, sigma0_linear


def calculate_rcs_calibrated(target_window_complex, noise_power, K):
    """
    Расчет калиброванной ЭПР для точечной цели
    """
    # Суммарная энергия цели в окне
    target_energy = np.sum(np.abs(target_window_complex) ** 2)

    # Вычитание энергии шума (шумовая мощность * количество пикселей)
    noise_energy = noise_power * target_window_complex.size
    calibrated_energy = target_energy - noise_energy

    # Расчет ЭПР
    rcs_linear = calibrated_energy / K
    rcs_db = 10 * np.log10(rcs_linear) if rcs_linear > 0 else -80

    return rcs_db, rcs_linear


def find_shadow_regions(radar_image_complex, min_region_size=50):
    """
    Поиск участков радиолокационной тени для оценки NESZ
    """
    radar_image_db = 20 * np.log10(np.abs(radar_image_complex) + 1e-12)

    # Порог для тени (значения ниже этого считаются тенями)
    shadow_threshold = np.percentile(radar_image_db, 5)

    shadow_regions = []
    h, w = radar_image_db.shape

    # Поиск связных областей с низкой яркостью
    shadow_mask = radar_image_db < shadow_threshold
    labeled_mask, num_features = ndimage.label(shadow_mask)

    for i in range(1, num_features + 1):
        region_mask = labeled_mask == i
        region_size = np.sum(region_mask)

        if region_size >= min_region_size:
            # Координаты региона
            region_coords = np.where(region_mask)
            y_min, y_max = np.min(region_coords[0]), np.max(region_coords[0])
            x_min, x_max = np.min(region_coords[1]), np.max(region_coords[1])

            # Средняя яркость региона
            region_brightness = np.mean(radar_image_db[region_mask])

            shadow_regions.append({
                'coords': (y_min, y_max, x_min, x_max),
                'size': region_size,
                'mean_brightness_db': region_brightness
            })

    return shadow_regions


СтерЖенёк(Ферритовый), [28.11.2025 2: 53]

def main():
    # Существующий код загрузки данных
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    detected_peaks = find_targets(radar_image_complex, MIN_DISTANCE, THRESHOLD_OFFSET_DB)
    T_synth, F_r_discr, Full_velocity, output_samp_rate, output_prf = read_radar_params_from_json(JSON_FILE_PATH)
    Noise_power = calculate_noise_threshold(radar_image_complex)

    # === ИЗМЕНЕНИЕ: Расчет шагов дискретизации ===
    d_r = c / (2 * output_samp_rate)  # Шаг по дальности
    d_az = Full_velocity / output_prf  # Шаг по азимуту

    # === ИЗМЕНЕНИЕ: Параметры калибровки ===
    # Для уголкового отражателя с ребром 0.5 м на частоте 9.6 ГГц
    sigma_ref = 268.5  # м² (расчетная ЭПР)

    # Координаты калибровочных целей (нужно задать или определить автоматически)
    # ВАЖНО: Замените на реальные координаты ваших калибровочных целей
    calibration_targets_coords = [
        (350, 1000),  # Пример координат калибровочной цели
        # Добавьте другие калибровочные цели при наличии
    ]

    # === ИЗМЕНЕНИЕ: Расчет калибровочного коэффициента ===
    K, E_mean_db = calculate_calibration_coefficient(
        radar_image_complex, calibration_targets_coords, sigma_ref, d_az, d_r
    )

    # === ИЗМЕНЕНИЕ: Расчет радиометрической чувствительности ===
    nesz_db, nesz_linear = calculate_nesz(Noise_power, K, d_az, d_r)

    # === ИЗМЕНЕНИЕ: Поиск теневых регионов для оценки NESZ ===
    shadow_regions = find_shadow_regions(radar_image_complex)

    targets_data = []

    for i, target_yx in enumerate(detected_peaks):
        # Выделение окна
        window = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        snr_db = calculate_snr_for_target(window, Noise_power)

        # === ИЗМЕНЕНИЕ: Расчет калиброванной ЭПР ===
        rcs_db, rcs_linear = calculate_rcs_calibrated(window, Noise_power, K)

        # === ИЗМЕНЕНИЕ: Расчет калиброванной УЭПР для центрального пикселя ===
        center_y, center_x = window.shape[0] // 2, window.shape[1] // 2
        center_power = np.abs(window[center_y, center_x]) ** 2
        sigma0_db, sigma0_linear = calculate_sigma0_calibrated(
            center_power, Noise_power, K, d_az, d_r
        )

        # Извлечение сечений (существующий код)
        horizontal_section, vertical_section = extract_sections(window)

        # Анализ сечений (существующий код)
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)

        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)

        # Формируем данные для отчета
        target_data = {
            'window_linear': window,
            'snr_db': snr_db,
            # === ИЗМЕНЕНИЕ: Добавляем калиброванные параметры ===
            'rcs_db': rcs_db,
            'rcs_linear': rcs_linear,
            'sigma0_db': sigma0_db,
            'sigma0_linear': sigma0_linear,
            # Существующие данные...
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

            СтерЖенёк(Ферритовый), [28.11.2025 2:53]
        'v_pslr': v_results.get('classical_pslr', -80),
        'v_i_pslr': v_results.get('integral_pslr', -80),
        'v_sinc_interp': v_results.get('sinc_interp'),
        'v_t_interp': v_results.get('t_interp'),
        }

        targets_data.append(target_data)

    # === ИЗМЕНЕНИЕ: Формирование отчета с радиометрическими параметрами ===
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Название голограммы: {HOLOGRAM_NAME}\n")
        f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
        f.write(f"Количество целей: {len(detected_peaks)}\n")
        f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
        f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
        f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")

        # Радиометрические параметры
        f.write("\n=== РАДИОМЕТРИЧЕСКИЕ ПАРАМЕТРЫ ===\n")
        f.write(f"Шаг по азимуту (d_az): {d_az:.6f} м\n")
        f.write(f"Шаг по дальности (d_r): {d_r:.6f} м\n")
        f.write(f"Расчетная ЭПР эталона (σ_ref): {sigma_ref:.2f} м²\n")
        f.write(f"Средняя энергия эталонов (E_mean): {E_mean_db:.2f} дБ\n")
        f.write(f"Калибровочный коэффициент (K): {K:.6f}\n")
        f.write(f"Радиометрическая чувствительность (NESZ): {nesz_db:.2f} дБ\n")
        f.write(f"Мощность шума: {Noise_power:.6f}\n")

        # Информация о теневых регионах
        f.write(f"\n=== ТЕНЕВЫЕ РЕГИОНЫ ДЛЯ ОЦЕНКИ NESZ ===\n")
        f.write(f"Количество теневых регионов: {len(shadow_regions)}\n")
        for i, region in enumerate(shadow_regions[:5]):  # Показываем первые 5
            f.write(f"Регион {i + 1}: размер={region['size']} пикс., УЭПР={region['mean_brightness_db']:.2f} дБ\n")

    # Существующий код визуализации и генерации отчета...
    # [остальная часть кода без изменений]