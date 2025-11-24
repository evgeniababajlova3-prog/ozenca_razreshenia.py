# === ПРОСТАЯ И ПОНЯТНАЯ ФУНКЦИЯ ОБНАРУЖЕНИЯ ЦЕЛЕЙ ===
def find_targets_by_snr(radar_image_complex, noise_power, min_distance=MIN_DISTANCE, min_snr_db=MIN_SNR_DB):
    """
    Простое обнаружение целей: находим пики, выделяем окно, считаем SNR в окне
    
    Parameters:
    -----------
    radar_image_complex : ndarray
        Комплексное радиолокационное изображение
    noise_power : float
        Мощность шума
    min_distance : int
        Минимальное расстояние между целями
    min_snr_db : float
        Минимальный SNR для обнаружения цели в дБ
    
    Returns:
    --------
    peaks_coords : list of tuples
        Координаты обнаруженных целей (y, x)
    snr_values : list
        Значения SNR для каждой обнаруженной цели
    """
    
    # 1. Находим все локальные максимумы амплитуды
    amplitude_image = np.abs(radar_image_complex)
    local_max = ndimage.maximum_filter(amplitude_image, size=min_distance) == amplitude_image
    candidate_peaks = np.where(local_max)
    candidate_coords = list(zip(candidate_peaks[0], candidate_peaks[1]))
    
    print(f"Найдено {len(candidate_coords)} кандидатов в цели")
    
    # 2. Фильтруем кандидатов по SNR в окне
    peaks_coords = []
    snr_values = []
    
    for y, x in candidate_coords:
        # Выделяем окно вокруг кандидата
        window_complex = extract_target_window(radar_image_complex, (y, x), WINDOW_SIZE)
        
        # Считаем средний SNR в окне (точно так же как при анализе цели)
        target_amplitudes = np.abs(window_complex).flatten()
        signal_power = np.mean(target_amplitudes ** 2)  # СРЕДНЯЯ мощность в окне
        snr_linear = signal_power / noise_power
        snr_db = 10 * np.log10(snr_linear)
        
        # Если SNR превышает порог - это цель
        if snr_db >= min_snr_db:
            peaks_coords.append((y, x))
            snr_values.append(snr_db)
    
    print(f"Обнаружено {len(peaks_coords)} целей с SNR > {min_snr_db} дБ")
    
    # Создаем простую карту SNR для визуализации (опционально)
    snr_map = np.zeros_like(amplitude_image)
    for (y, x), snr in zip(peaks_coords, snr_values):
        snr_map[y, x] = snr
    
    return peaks_coords, snr_values, snr_map

# === ОБНОВЛЕННАЯ ФУНКЦИЯ calculate_snr_for_target (теперь простая) ===
def calculate_snr_for_target(window_complex, noise_power):
    """
    Простой расчет SNR для цели - тот же алгоритм, что при обнаружении
    
    Parameters:
    -----------
    window_complex : ndarray
        Комплексное окно вокруг цели  
    noise_power : float
        Мощность шума
        
    Returns:
    --------
    snr_db : float
        SNR в дБ
    """
    target_amplitudes = np.abs(window_complex).flatten()
    signal_power = np.mean(target_amplitudes ** 2)
    snr_linear = signal_power / noise_power
    snr_db = 10 * np.log10(snr_linear)
    return snr_db

# === УПРОЩЕННАЯ ИНТЕГРАЦИЯ В MAIN() ===
def main():
    # ... загрузка данных ...
    
    # Расчет параметров шума
    noise_rms, noise_power, noise_amplitudes = calculate_noise_threshold(radar_image_complex)
    
    # Простое обнаружение целей
    print("=== ПРОСТОЕ ОБНАРУЖЕНИЕ ЦЕЛЕЙ ПО SNR В ОКНЕ ===")
    detected_peaks, detected_snrs, snr_map = find_targets_by_snr(
        radar_image_complex, noise_power, MIN_DISTANCE, MIN_SNR_DB
    )
    
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)
    
    targets_data = []
    
    for i, (target_yx, snr_db) in enumerate(zip(detected_peaks, detected_snrs)):
        # Выделяем окно для цели
        window_complex = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        window_amplitude = np.abs(window_complex)
        
        # SNR уже рассчитан при обнаружении, используем его
        # Можно пересчитать для проверки:
        snr_db_check = calculate_snr_for_target(window_complex, noise_power)
        
        print(f"Цель {i+1}: SNR={snr_db:.2f} дБ (проверка: {snr_db_check:.2f} дБ)")
        
        # Извлечение сечений
        horizontal_section, vertical_section = extract_sections(window_amplitude)
        
        # Анализ сечений
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)
        
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)
        
        # Формируем данные для отчета
        target_data = {
            'window_linear': window_amplitude,
            'window_complex': window_complex,
            'snr_db': snr_db,  # Используем SNR, рассчитанный при обнаружении
            # ... остальные поля без изменений ...
        }
        targets_data.append(target_data)
    
    # ... остальная часть main ...    local_max_snr = ndimage.maximum_filter(snr_map, size=min_distance) == snr_map
    
    # Применяем порог по SNR
    above_snr_threshold = snr_map > min_snr_db
    
    # Комбинируем условия: локальные максимумы SNR с достаточным SNR
    detected = local_max_snr & above_snr_threshold
    
    # Получаем координаты целей
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))
    
    print(f"Обнаружено {len(peaks_coords)} целей с усредненным SNR > {min_snr_db} дБ")
    print(f"Размер окна усреднения мощности: {power_window_size}")
    print(f"Мощность шума: {noise_power:.6f}")
    
    # Отладочная информация
    if len(peaks_coords) > 0:
        snr_values = [snr_map[y, x] for y, x in peaks_coords]
        power_values = [averaged_power[y, x] for y, x in peaks_coords]
        
        print(f"SNR целей: min={min(snr_values):.2f} дБ, max={max(snr_values):.2f} дБ, mean={np.mean(snr_values):.2f} дБ")
        print(f"Мощность целей: min={min(power_values):.6f}, max={max(power_values):.6f}")
    
    return peaks_coords, snr_map

# === ВСПОМОГАТЕЛЬНАЯ ФУНКЦИЯ ДЛЯ РАСЧЕТА SNR ЦЕЛИ ===
def calculate_snr_for_target(window_complex, noise_power):
    """
    Расчет среднего SNR для цели (согласовано с обнаружением)
    
    Parameters:
    -----------
    window_complex : ndarray
        Комплексное окно вокруг цели  
    noise_power : float
        СРЕДНЯЯ мощность шума
        
    Returns:
    --------
    snr_db : float
        SNR в дБ
    """
    
    # СРЕДНЯЯ мощность сигнала в окне цели
    target_amplitudes = np.abs(window_complex).flatten()
    signal_power = np.mean(target_amplitudes ** 2)  # СРЕДНЯЯ мощность!
    
    # SNR в линейной области
    snr_linear = signal_power / noise_power
    
    # SNR в дБ
    snr_db = 10 * np.log10(snr_linear)
    
    return snr_db

# === ОБНОВЛЕННЫЙ ВЫЗОВ В MAIN() ===
def main():
    # ... загрузка данных ...
    
    # Расчет параметров шума (должна возвращать СРЕДНЮЮ мощность шума)
    noise_rms, noise_power, noise_amplitudes = calculate_noise_threshold(radar_image_complex)
    
    # === ИСПРАВЛЕННЫЙ ВЫЗОВ: Обнаружение по усредненной мощности ===
    detected_peaks, snr_map = find_targets_by_snr(
        radar_image_complex, 
        noise_power, 
        MIN_DISTANCE, 
        MIN_SNR_DB, 
        power_window_size=(5, 5)  # Можно настроить: (3,3), (7,7), и т.д.
    )
    
    # Обработка целей...
    for i, target_yx in enumerate(detected_peaks):
        window_complex = extract_target_window(radar_image    
    # # === ПРАВИЛЬНАЯ ФУНКЦИЯ ОБНАРУЖЕНИЯ ЦЕЛЕЙ ПО SNR ===
def find_targets_by_snr(radar_image_complex, noise_power, min_distance=MIN_DISTANCE, min_snr_db=MIN_SNR_DB):
    """
    Обнаружение целей на основе локальных максимумов карты SNR
    
    Parameters:
    -----------
    radar_image_complex : ndarray
        Комплексное радиолокационное изображение
    noise_power : float
        Мощность шума
    min_distance : int
        Минимальное расстояние между целями
    min_snr_db : float
        Минимальный SNR для обнаружения цели в дБ
    
    Returns:
    --------
    peaks_coords : list of tuples
        Координаты обнаруженных целей (y, x)
    snr_map : ndarray
        Карта SNR для всего изображения
    """
    
    # Создаем карту амплитуд
    amplitude_image = np.abs(radar_image_complex)
    
    # === ИСПРАВЛЕНИЕ: Правильный расчет карты SNR ===
    # SNR = 10 * log10(мощность_сигнала / мощность_шума)
    power_image = amplitude_image ** 2
    snr_linear = power_image / noise_power
    snr_map = 10 * np.log10(np.maximum(snr_linear, 1e-12))  # защита от log(0)
    
    # === ИСПРАВЛЕНИЕ: Ищем локальные максимумы НА КАРТЕ SNR, а не амплитуды ===
    # Применяем фильтр максимумов к карте SNR
    local_max_snr = ndimage.maximum_filter(snr_map, size=min_distance) == snr_map
    
    # Применяем порог по SNR
    above_snr_threshold = snr_map > min_snr_db
    
    # Комбинируем условия: локальные максимумы SNR с достаточным SNR
    detected = local_max_snr & above_snr_threshold
    
    # Получаем координаты целей
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))
    
    print(f"Обнаружено {len(peaks_coords)} целей с SNR > {min_snr_db} дБ")
    print(f"Мощность шума: {noise_power:.6f}")
    
    # Отладочная информация
    if len(peaks_coords) > 0:
        snr_values = [snr_map[y, x] for y, x in peaks_coords]
        print(f"SNR целей: min={min(snr_values):.2f} дБ, max={max(snr_values):.2f} дБ, mean={np.mean(snr_values):.2f} дБ")
    
    return peaks_coords, snr_map

# === АЛЬТЕРНАТИВНАЯ ФУНКЦИЯ: Обнаружение по среднему SNR в окне ===
def find_targets_by_mean_snr(radar_image_complex, noise_power, window_size=WINDOW_SIZE, min_distance=MIN_DISTANCE, min_snr_db=MIN_SNR_DB):
    """
    Обнаружение целей на основе среднего SNR в окне вокруг локальных максимумов амплитуды
    
    Parameters:
    -----------
    radar_image_complex : ndarray
        Комплексное радиолокационное изображение
    noise_power : float
        Мощность шума
    window_size : tuple
        Размер окна для расчета среднего SNR
    min_distance : int
        Минимальное расстояние между целями
    min_snr_db : float
        Минимальный SNR для обнаружения цели в дБ
    
    Returns:
    --------
    peaks_coords : list of tuples
        Координаты обнаруженных целей (y, x)
    snr_map : ndarray
        Карта SNR для всего изображения
    """
    
    amplitude_image = np.abs(radar_image_complex)
    
    # Создаем карту SNR (как в предыдущей функции)
    power_image = amplitude_image ** 2
    snr_linear = power_image / noise_power
    snr_map = 10 * np.log10(np.maximum(snr_linear, 1e-12))
    
    # Находим локальные максимумы амплитуды
    local_max_amplitude = ndimage.maximum_filter(amplitude_image, size=min_distance) == amplitude_image
    
    peaks = np.where(local_max_amplitude)
    candidate_coords = list(zip(peaks[0], peaks[1]))
    
    peaks_coords = []
    
    for y, x in candidate_coords:
        # Выделяем окно вокруг кандидата
        window_complex = extract_target_window(radar_image_complex, (y, x), window_size)
        
        # Рассчитываем средний SNR в окне (точно так же как в calculate_snr_for_target)
        target_amplitudes = np.abs(window_complex).flatten()
        signal_power = np.mean(target_amplitudes ** 2)
        snr_linear = signal_power / noise_power
        snr_db = 10 * np.log10(snr_linear)
        
        # Если средний SNR превышает порог - это цель
        if snr_db >= min_snr_db:
            peaks_coords.append((y, x))
    
    print(f"Обнаружено {len(peaks_coords)} целей со средним SNR в окне > {min_snr_db} дБ")
    
    return peaks_coords, snr_map

# === ИСПРАВЛЕННАЯ ИНТЕГРАЦИЯ В ОСНОВНУЮ ФУНКЦИЮ main() ===

def main():
    # Загрузка данных
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    
    # Расчет параметров шума
    print("=== РАСЧЕТ ПАРАМЕТРОВ ШУМА ДЛЯ SNR ===")
    noise_rms, noise_power, noise_amplitudes = calculate_noise_threshold(radar_image_complex)
    
    # === ВЫБОР МЕТОДА ОБНАРУЖЕНИЯ ===
    # Вариант 1: Обнаружение по локальным максимумам SNR (рекомендуется)
    print("=== ОБНАРУЖЕНИЕ ЦЕЛЕЙ ПО ЛОКАЛЬНЫМ МАКСИМУМАМ SNR ===")
    detected_peaks, snr_map = find_targets_by_snr(radar_image_complex, noise_power, MIN_DISTANCE, MIN_SNR_DB)
    
    # Вариант 2: Обнаружение по среднему SNR в окне (раскомментируйте если нужно)
    # print("=== ОБНАРУЖЕНИЕ ЦЕЛЕЙ ПО СРЕДНЕМУ SNR В ОКНЕ ===")
    # detected_peaks, snr_map = find_targets_by_mean_snr(radar_image_complex, noise_power, WINDOW_SIZE, MIN_DISTANCE, MIN_SNR_DB)
    
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)
    
    targets_data = []
    
    for i, target_yx in enumerate(detected_peaks):
        # Выделяем комплексное окно для цели
        window_complex = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        window_amplitude = np.abs(window_complex)
        
        # === РАСЧЕТ SNR ДЛЯ ЦЕЛИ (ТОЧНО ТАК ЖЕ КАК ПРИ ОБНАРУЖЕНИИ) ===
        if 'find_targets_by_mean_snr' in globals() and globals()['find_targets_by_mean_snr'].__code__ == find_targets_by_mean_snr.__code__:
            # Если использовали обнаружение по среднему SNR в окне
            target_amplitudes = np.abs(window_complex).flatten()
            signal_power = np.mean(target_amplitudes ** 2)
            snr_linear = signal_power / noise_power
            snr_db = 10 * np.log10(snr_linear)
        else:
            # Если использовали обнаружение по локальным максимумам SNR
            # Берем SNR из карты в точке максимума (это тот SNR по которому обнаружили цель)
            snr_db = snr_map[target_yx[0], target_yx[1]]
        
        # Дополнительно: пиковый SNR (для информации)
        peak_amplitude = np.max(window_amplitude)
        snr_peak_db = 20 * np.log10(peak_amplitude / noise_rms)
        
        print(f"Цель {i+1}: SNR={snr_db:.2f} дБ, SNR_пиковый={snr_peak_db:.2f} дБ")
        
        # Извлечение сечений из амплитудного окна
        horizontal_section, vertical_section = extract_sections(window_amplitude)
        
        # Анализ сечений
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)
        
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)
        
        # Формируем данные для отчета
        target_data = {
            'window_linear': window_amplitude,
            'window_complex': window_complex,
            'snr_db': snr_db,                    # SNR по которому цель была обнаружена
            'snr_peak_db': snr_peak_db,          # Пиковый SNR (дополнительно)
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
    
    # ... остальная часть main без изменений ... порог по SNR
    above_snr_threshold = snr_map > min_snr_db
    
    # Комбинируем условия: локальные максимумы с достаточным SNR
    detected = local_max & above_snr_threshold
    
    # Получаем координаты целей
    peaks = np.where(detected)
    peaks_coords = list(zip(peaks[0], peaks[1]))
    
    print(f"Обнаружено {len(peaks_coords)} целей с SNR > {min_snr_db} дБ")
    print(f"Мощность шума: {noise_power:.6f}, RMS шума: {noise_rms:.6f}")
    
    # Отладочная информация
    if len(peaks_coords) > 0:
        snr_values = [snr_map[y, x] for y, x in peaks_coords]
        print(f"SNR целей: min={min(snr_values):.2f} дБ, max={max(snr_values):.2f} дБ, mean={np.mean(snr_values):.2f} дБ")
    
    return peaks_coords, snr_map

def find_targets_by_peak_snr(radar_image_complex, noise_rms, min_distance=MIN_DISTANCE, min_snr_db=MIN_SNR_DB):
    """
    Альтернативный метод обнаружения по пиковому SNR
    (отношение максимальной амплитуды в окне к RMS шума)
    """
    
    amplitude_image = np.abs(radar_image_complex)
    
    # Находим локальные максимумы
    local_max = ndimage.maximum_filter(amplitude_image, size=min_distance) == amplitude_image
    
    peaks = np.where(local_max)
    peaks_coords = []
    
    for y, x in zip(peaks[0], peaks[1]):
        # Вычисляем пиковый SNR для кандидата
        peak_amplitude = amplitude_image[y, x]
        peak_snr_db = 20 * np.log10(peak_amplitude / noise_rms)
        
        if peak_snr_db >= min_snr_db:
            peaks_coords.append((y, x))
    
    print(f"Обнаружено {len(peaks_coords)} целей с пиковым SNR > {min_snr_db} дБ")
    
    return peaks_coords

# === ОБНОВЛЕННАЯ ФУНКЦИЯ calculate_noise_threshold ===
def calculate_noise_threshold(radar_image_complex, snr_threshold_db=MIN_SNR_DB):
    """
    Расчет параметров шума и порога обнаружения на основе SNR
    
    Returns:
    --------
    noise_rms : float
        RMS шума
    noise_power : float
        Мощность шума  
    noise_amplitudes : ndarray
        Амплитуды шумовой области
    """
    
    # Находим шумовую область
    noise_window, noise_amplitudes, noise_phases, match_quality = find_rayleigh_uniform_region(radar_image_complex)
    
    if noise_amplitudes is None:
        raise ValueError("Не удалось найти шумовую область для расчета SNR")
    
    # Вычисляем параметры шума
    noise_power = np.mean(noise_amplitudes ** 2)
    noise_rms = np.sqrt(noise_power)
    
    print(f"Параметры шума: мощность={noise_power:.6f}, RMS={noise_rms:.6f}")
    print(f"Порог обнаружения: SNR > {snr_threshold_db} дБ")
    
    return noise_rms, noise_power, noise_amplitudes

# === ИНТЕГРАЦИЯ В ОСНОВНУЮ ФУНКЦИЮ main() ===

def main():
    # Загрузка данных
    radar_image_complex, radar_image = load_radar_image_from_hdf5(HDF5_FILE_PATH)
    
    # === ИЗМЕНЕНИЕ: Расчет параметров шума для SNR ===
    print("=== РАСЧЕТ ПАРАМЕТРОВ ШУМА ДЛЯ SNR ===")
    noise_rms, noise_power, noise_amplitudes = calculate_noise_threshold(radar_image_complex, MIN_SNR_DB)
    
    # === ИЗМЕНЕНИЕ: Обнаружение целей по SNR ===
    print("=== ОБНАРУЖЕНИЕ ЦЕЛЕЙ ПО SNR ===")
    detected_peaks, snr_map = find_targets_by_snr(radar_image_complex, noise_amplitudes, MIN_DISTANCE, MIN_SNR_DB)
    
    # Альтернативный метод (раскомментируйте если нужно):
    # detected_peaks = find_targets_by_peak_snr(radar_image_complex, noise_rms, MIN_DISTANCE, MIN_SNR_DB)
    
    T_synth, F_r_discr, Full_velocity = read_radar_params_from_json(JSON_FILE_PATH)
    
    targets_data = []
    
    for i, target_yx in enumerate(detected_peaks):
        # Выделяем комплексное окно для цели
        window_complex = extract_target_window(radar_image_complex, target_yx, WINDOW_SIZE)
        window_amplitude = np.abs(window_complex)
        
        # === ДОБАВЛЕНИЕ: Расчет различных метрик SNR ===
        # Средний SNR по окну
        window_power = np.mean(window_amplitude ** 2)
        snr_db = 10 * np.log10(window_power / noise_power)
        
        # Пиковый SNR
        peak_amplitude = np.max(window_amplitude)
        snr_peak_db = 20 * np.log10(peak_amplitude / noise_rms)
        
        # SNR в точке максимума из карты SNR
        snr_at_peak = snr_map[target_yx[0], target_yx[1]]
        
        print(f"Цель {i+1}: SNR_окно={snr_db:.2f} дБ, SNR_пик={snr_peak_db:.2f} дБ, SNR_карта={snr_at_peak:.2f} дБ")
        
        # Извлечение сечений из амплитудного окна
        horizontal_section, vertical_section = extract_sections(window_amplitude)
        
        # Анализ сечений
        t_h, h_signal_db, h_signal_linear = generate_sinc_signal_from_section(horizontal_section, WINDOW_SIZE[0])
        h_results = analiz_sechenia(t_h, h_signal_db, WINDOW_SIZE[0] / 2)
        
        t_v, v_signal_db, v_signal_linear = generate_sinc_signal_from_section(vertical_section, WINDOW_SIZE[1])
        v_results = analiz_sechenia(t_v, v_signal_db, WINDOW_SIZE[1] / 2)
        
        # Формируем данные для отчета
        target_data = {
            'window_linear': window_amplitude,
            'window_complex': window_complex,
            'snr_db': snr_db,                    # Средний SNR по окну
            'snr_peak_db': snr_peak_db,          # Пиковый SNR
            'snr_map_db': snr_at_peak,           # SNR из карты в точке максимума
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
    
    # === ДОПОЛНИТЕЛЬНАЯ ВИЗУАЛИЗАЦИЯ: Карта SNR ===
    plt.figure(figsize=(12, 8))
    plt.imshow(snr_map, cmap="jet", aspect='auto')
    plt.colorbar(label='SNR (дБ)')
    plt.title('Карта SNR радиолокационного изображения')
    
    # Отмечаем обнаруженные цели
    for i, (y, x) in enumerate(detected_peaks):
        plt.plot(x, y, 's', markersize=12, markeredgewidth=2, 
                markeredgecolor='red', markerfacecolor='none', linestyle='none',
                label='Цели' if i == 0 else "")
    
    plt.legend()
    snr_map_path = os.path.join(result_folder, "snr_map.png")
    plt.savefig(snr_map_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # ... остальная часть main без изменений до генерации отчетов ...
    
    # === ИЗМЕНЕНИЕ В ОТЧЕТЕ: Добавляем информацию о SNR и методе обнаружения ===
    
    # В текстовом отчете о параметрах
    with open(os.path.join(result_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Название голограммы: {HOLOGRAM_NAME}\n")
        f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
        f.write(f"Количество целей: {len(detected_peaks)}\n")
        f.write(f"Метод обнаружения: по SNR (порог: {MIN_SNR_DB} дБ)\n")
        f.write(f"Мощность шума: {noise_power:.6f}\n")
        f.write(f"RMS шума: {noise_rms:.6f}\n")
        f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
        f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
        f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")
    
    # В отчетах для каждой цели
    for i, target_data in enumerate(targets_data):
        target_id = i + 1
        target_coords = detected_peaks[i]
        
        with  right_minima = minima_indices[minima_indices > main_peak_idx]

    # ИЗМЕНЕНИЕ 1: Определяем границы главного лепестка с учетом краев
    if len(left_minima) == 0:
        # Если нет левых минимумов, берем начало сигнала как левую границу
        left_zero_idx = 0
        print("Левый ноль не найден, используем начало сигнала")
    else:
        left_zero_idx = left_minima[-1]  # последний ноль слева

    if len(right_minima) == 0:
        # Если нет правых минимумов, берем конец сигнала как правую границу
        right_zero_idx = len(signal_linear) - 1
        print("Правый ноль не найден, используем конец сигнала")
    else:
        right_zero_idx = right_minima[0]  # первый ноль справа

    # Главный лепесток - между этими нулями
    main_lobe_indices = range(left_zero_idx, right_zero_idx + 1)
    main_lobe = signal_linear[main_lobe_indices]
    t_main = t[main_lobe_indices]

    # ИЗМЕНЕНИЕ 2: Боковые лепестки - все что осталось, даже если только с одной стороны
    sidelobe_indices = []
    
    # Левые боковые лепестки (если есть)
    if left_zero_idx > 0:
        sidelobe_indices.extend(np.arange(0, left_zero_idx))
    
    # Правые боковые лепестки (если есть)  
    if right_zero_idx < len(signal_linear) - 1:
        sidelobe_indices.extend(np.arange(right_zero_idx + 1, len(signal_linear)))
    
    sidelobe_indices = np.array(sidelobe_indices)
    
    if len(sidelobe_indices) == 0:
        print("Боковые лепестки не обнаружены (возможно, весь сигнал - главный лепесток)")
        return -80, -80

    sidelobes = signal_linear[sidelobe_indices]
    t_sidelobes = t[sidelobe_indices]

    # ИЗМЕНЕНИЕ 3: Классический УБЛ - отношение максимального бокового лепестка к главному
    if len(sidelobes) > 0:
        max_sidelobe = np.max(sidelobes)
        classical_pslr_db = 20 * np.log10(max_sidelobe / main_peak_val)
    else:
        classical_pslr_db = -80

    # ИЗМЕНЕНИЕ 4: Интегральный УБЛ - отношение мощностей (работает даже с частичными данными)
    if len(main_lobe) > 0 and len(sidelobes) > 0:
        power_main = trapezoid(main_lobe ** 2, t_main)
        power_sidelobes = trapezoid(sidelobes ** 2, t_sidelobes)
        
        # Если один из боковых лепестков отсутствует, делаем поправку
        left_sidelobes_present = left_zero_idx > 0
        right_sidelobes_present = right_zero_idx < len(signal_linear) - 1
        
        if not left_sidelobes_present or not right_sidelobes_present:
            # Если отсутствует один из боковых лепестков, предполагаем симметрию для оценки
            print(f"УБЛ рассчитан для частичных данных: левый={left_sidelobes_present}, правый={right_sidelobes_present}")
            # Можно добавить коэффициент коррекции, но пока просто используем как есть
        
        integral_pslr_db = 10 * np.log10(power_sidelobes / power_main)
    else:
        integral_pslr_db = -80

    # ИЗМЕНЕНИЕ 5: Дополнительная диагностика
    print(f"Главный лепесток: отсчеты {left_zero_idx}-{right_zero_idx}, "
          f"боковые лепестки: {len(sidelobes)} отсчетов")
    
    # Визуализация для отладки (можно закомментировать в финальной версии)
    if False:  # Поставьте True для отладкиimport numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import firwin, lfilter
import warnings
warnings.filterwarnings('ignore')

def sinc_interpolation_wong(t_original, signal_db, interp_factor=4):
    """Sinc-интерполяция по алгоритму Wong 2005"""
    
    kernel_size = 8
    subsamples = 16
    
    def generate_coefficients(kernel_size, subsamples):
        coefficients = np.zeros((subsamples, kernel_size))
        window = np.kaiser(kernel_size, beta=2.5)

        for phase in range(subsamples):
            shift = phase / subsamples
            for k in range(kernel_size):
                pos = (k - (kernel_size-1)//2) - shift
                if abs(pos) < 1e-12:
                    sinc_val = 1.0
                else:
                    sinc_val = np.sin(np.pi * pos) / (np.pi * pos)
                coefficients[phase, k] = sinc_val * window[k]

            sum_coeff = np.sum(coefficients[phase])
            if abs(sum_coeff) > 1e-12:
                coefficients[phase] /= sum_coeff
        return coefficients

    # Основная функция интерполяции
    coefficients_table = generate_coefficients(kernel_size, subsamples)
    signal_linear = 10 ** (signal_db / 20)
    dt = t_original[1] - t_original[0]
    
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    signal_interp = np.zeros_like(t_interp)
    half_kernel = kernel_size // 2

    for i, t_point in enumerate(t_interp):
        idx_center = int(np.floor((t_point - t_original[0]) / dt))
        fractional = (t_point - t_original[idx_center]) / dt
        
        phase_idx = int(fractional * subsamples)
        phase_idx = max(0, min(phase_idx, subsamples - 1))
        coeffs = coefficients_table[phase_idx]

        start_idx = idx_center - half_kernel
        end_idx = start_idx + kernel_size

        if start_idx < 0 or end_idx > len(signal_linear):
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
        else:
            points = signal_linear[start_idx:end_idx]
            signal_interp[i] = np.sum(points * coeffs)

    signal_interp_db = 20 * np.log10(np.abs(signal_interp) + 1e-12)
    return t_interp, signal_interp_db

def fourier_interpolation_basic(t_original, signal_db, interp_factor=4):
    """Базовая Фурье-интерполяция с дополнением нулями"""
    
    signal_linear = 10 ** (signal_db / 20)
    N_original = len(signal_linear)
    
    spectrum = np.fft.fft(signal_linear)
    N_new = N_original * interp_factor
    new_spectrum = np.zeros(N_new, dtype=complex)
    
    half_original = N_original // 2
    new_spectrum[:half_original] = spectrum[:half_original]
    new_spectrum[-half_original:] = spectrum[-half_original:]
    new_spectrum *= interp_factor
    
    signal_interp_linear = np.fft.ifft(new_spectrum).real
    t_interp = np.linspace(t_original[0], t_original[-1], N_new)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def fourier_interpolation_hann(t_original, signal_db, interp_factor=4):
    """Фурье-интерполяция с оконной функцией Hann"""
    
    signal_linear = 10 ** (signal_db / 20)
    N_original = len(signal_linear)
    
    # Оконная функция
    window = np.hanning(N_original)
    signal_windowed = signal_linear * window
    
    spectrum = np.fft.fft(signal_windowed)
    N_new = N_original * interp_factor
    new_spectrum = np.zeros(N_new, dtype=complex)
    
    half_original = N_original // 2
    new_spectrum[:half_original] = spectrum[:half_original]
    new_spectrum[-half_original:] = spectrum[-half_original:]
    new_spectrum *= interp_factor
    
    signal_interp_linear = np.fft.ifft(new_spectrum).real
    
    # Компенсация оконной функции
    window_interp = np.interp(
        np.linspace(0, N_original-1, N_new),
        np.arange(N_original),
        window
    )
    
    signal_interp_linear = np.divide(
        signal_interp_linear, 
        window_interp + 1e-12,
        out=signal_interp_linear,
        where=window_interp > 1e-6
    )
    
    t_interp = np.linspace(t_original[0], t_original[-1], N_new)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def cubic_spline_interpolation(t_original, signal_db, interp_factor=4):
    """Кубическая сплайн-интерполяция"""
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    
    cs = CubicSpline(t_original, signal_linear)
    signal_interp_linear = cs(t_interp)
    
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    return t_interp, signal_interp_db

def linear_interpolation(t_original, signal_db, interp_factor=4):
    """Линейная интерполяция"""
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    
    signal_interp_linear = np.interp(t_interp, t_original, signal_linear)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def lagrange_interpolation(t_original, signal_db, interp_factor=4, order=3):
    """Полиномиальная интерполяция методом Лагранжа"""
    
    def lagrange_polynomial(x, x_points, y_points):
        total = 0
        n = len(x_points)
        for i in range(n):
            xi, yi = x_points[i], y_points[i]
            product = 1
            for j in range(n):
                if i != j:
                    xj = x_points[j]
                    product *= (x - xj) / (xi - xj)
            total += yi * product
        return total
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    signal_interp_linear = np.zeros_like(t_interp)
    
    # Применяем скользящее окно Лагранжа
    for i, t_point in enumerate(t_interp):
        # Находим ближайшие точки
        distances = np.abs(t_original - t_point)
        nearest_indices = np.argsort(distances)[:order+1]
        nearest_t = t_original[nearest_indices]
        nearest_y = signal_linear[nearest_indices]
        
        signal_interp_linear[i] = lagrange_polynomial(t_point, nearest_t, nearest_y)
    
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    return t_interp, signal_interp_db

def fir_filter_interpolation(t_original, signal_db, interp_factor=4):
    """Интерполяция через FIR фильтр"""
    
    signal_linear = 10 ** (signal_db / 20)
    
    # Создаем FIR фильтр для интерполяции
    numtaps = 31
    cutoff = 0.8 / interp_factor
    fir_coeff = firwin(numtaps, cutoff)
    
    # Увеличиваем частоту дискретизации
    upsampled = np.zeros(len(signal_linear) * interp_factor)
    upsampled[::interp_factor] = signal_linear
    
    # Применяем FIR фильтр
    signal_interp_linear = lfilter(fir_coeff, 1.0, upsampled)
    
    # Обрезаем переходный процесс
    signal_interp_linear = signal_interp_linear[numtaps//2:]
    
    t_interp = np.linspace(t_original[0], t_original[-1], len(signal_interp_linear))
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def calculate_mse(t_original, signal_original_db, t_interp, signal_interp_db):
    """Вычисление среднеквадратичной ошибки между интерполированным и теоретическим сигналом"""
    
    # Создаем теоретический sinc сигнал на интерполированной сетке
    t_theoretical = np.linspace(t_original[0], t_original[-1], len(t_interp))
    ideal_sinc = np.sinc(t_theoretical)
    ideal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Вычисляем MSE в dB
    mse_db = np.mean((signal_interp_db - ideal_db) ** 2)
    
    # Также вычисляем MSE в линейной области для сравнения
    signal_interp_linear = 10 ** (signal_interp_db / 20)
    ideal_linear = np.abs(np.sinc(t_theoretical))
    mse_linear = np.mean((signal_interp_linear - ideal_linear) ** 2)
    
    return mse_db, mse_linear

def test_interpolation_methods():
    """Тестирование всех методов интерполяции на идеальном sinc-сигнале"""
    
    # Создаем идеальный sinc-сигнал для тестирования
    t_original = np.linspace(-3, 3, 30)  # 30 точек
    ideal_sinc = np.sinc(t_original)
    signal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Список методов для тестирования
    methods = [
        ('Sinc (Wong)', sinc_interpolation_wong),
        ('FFT (базовый)', fourier_interpolation_basic),
        ('FFT + Hann', fourier_interpolation_hann),
        ('Кубический сплайн', cubic_spline_interpolation),
        ('Линейная', linear_interpolation),
        ('Лагранж (порядок 3)', lambda t, s: lagrange_interpolation(t, s, 4, 3)),
        ('FIR фильтр', fir_filter_interpolation)
    ]
    
    # Результаты
    results = []
    
    # Создаем график для сравнения
    plt.figure(figsize=(15, 10))
    
    # Рисуем теоретический sinc
    t_dense = np.linspace(-3, 3, 1000)
    ideal_dense = np.sinc(t_dense)
    ideal_dense_db = 20 * np.log10(np.abs(ideal_dense) + 1e-12)
    plt.plot(t_dense, ideal_dense_db, 'k-', linewidth=3, label='Теоретический sinc', alpha=0.7)
    
    # Рисуем исходные точки
    plt.plot(t_original, signal_db, 'ko', markersize=8, label='Исходные точки', alpha=0.8)
    
    # Тестируем каждый метод
    colors = plt.cm.tab10(np.linspace(0, 1, len(methods)))
    
    for idx, (method_name, method_func) in enumerate(methods):
        try:
            # Применяем интерполяцию
            t_interp, signal_interp = method_func(t_original, signal_db, 4)
            
            # Вычисляем MSE
            mse_db, mse_linear = calculate_mse(t_original, signal_db, t_interp, signal_interp)
            
            # Сохраняем результаты
            results.append({
                'method': method_name,
                'mse_db': mse_db,
                'mse_linear': mse_linear,
                't_interp': t_interp,
                'signal_interp': signal_interp
            })
            
            # Рисуем интерполированный сигнал
            plt.plot(t_interp, signal_interp, 
                    color=colors[idx], 
                    linewidth=2, 
                    label=f'{method_name} (MSE: {mse_db:.4f})',
                    alpha=0.8)
            
            print(f"{method_name:<20} | MSE (dB): {mse_db:.6f} | MSE (linear): {mse_linear:.6f}")
            
        except Exception as e:
            print(f"{method_name:<20} | ОШИБКА: {str(e)}")
            results.append({
                'method': method_name,
                'mse_db': float('inf'),
                'mse_linear': float('inf'),
                'error': str(e)
            })
    
    # Настройка графика
    plt.xlabel('Время')
    plt.ylabel('Амплитуда (dB)')
    plt.title('Сравнение методов интерполяции на sinc-сигнале')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.ylim(-50, 5)
    plt.tight_layout()
    plt.show()
    
    # Сортируем методы по качеству (по MSE в dB)
    results_sorted = sorted([r for r in results if 'mse_db' in r], key=lambda x: x['mse_db'])
    
    # Выводим таблицу результатов
    print("\n" + "="*70)
    print("РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ МЕТОДОВ ИНТЕРПОЛЯЦИИ")
    print("="*70)
    print(f"{'Метод':<25} {'MSE (dB)':<15} {'MSE (linear)':<15}")
    print("-"*70)
    
    for result in results_sorted:
        print(f"{result['method']:<25} {result['mse_db']:<15.6f} {result['mse_linear']:<15.6f}")
    
    # Выводим лучший метод
    if results_sorted:
        best_method = results_sorted[0]
        print("-"*70)
        print(f"ЛУЧШИЙ МЕТОД: {best_method['method']}")
        print(f"MSE (dB): {best_method['mse_db']:.6f}")
        print(f"MSE (linear): {best_method['mse_linear']:.6f}")
    
    return results_sorted

def test_with_noisy_signal():
    """Дополнительный тест с зашумленным сигналом"""
    
    print("\n" + "="*70)
    print("ТЕСТ С ЗАШУМЛЕННЫМ SINC-СИГНАЛОМ")
    print("="*70)
    
    # Создаем зашумленный sinc-сигнал
    t_original = np.linspace(-3, 3, 30)
    ideal_sinc = np.sinc(t_original)
    
    # Добавляем шум
    np.random.seed(42)
    noise = np.random.normal(0, 0.05, len(ideal_sinc))
    noisy_sinc = ideal_sinc + noise
    signal_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
    
    # Тестируем только лучшие методы
    methods = [
        ('Sinc (Wong)', sinc_interpolation_wong),
        ('FFT + Hann', fourier_interpolation_hann),
        ('Кубический сплайн', cubic_spline_interpolation),
        ('FIR фильтр', fir_filter_interpolation)
    ]
    
    results = []
    
    for method_name, method_func in methods:
        try:
            t_interp, signal_interp = method_func(t_original, signal_db, 4)
            
            # Вычисляем MSE относительно идеального sinc (без шума)
            t_theoretical = np.linspace(t_original[0], t_original[-1], len(t_interp))
            ideal_sinc_dense = np.sinc(t_theoretical)
            ideal_db = 20 * np.log10(np.abs(ideal_sinc_dense) + 1e-12)
            
            mse_db = np.mean((signal_interp - ideal_db) ** 2)
            
            signal_interp_linear = 10 ** (signal_interp / 20)
            ideal_linear = np.abs(np.sinc(t_theoretical))
            mse_linear = np.mean((signal_interp_linear - ideal_linear) ** 2)
            
            results.append({
                'method': method_name,
                'mse_db': mse_db,
                'mse_linear': mse_linear
            })
            
            print(f"{method_name:<20} | MSE (dB): {mse_db:.6f} | MSE (linear): {mse_linear:.6f}")
            
        except Exception as e:
            print(f"{method_name:<20} | ОШИБКА: {str(e)}")
    
    # Сортируем по качеству
    results_sorted = sorted(results, key=lambda x: x['mse_db'])
    
    print("\nЛучший метод для зашумленного сигнала:")
    print(f"{results_sorted[0]['method']} с MSE (dB) = {results_sorted[0]['mse_db']:.6f}")
    
    return results_sorted

def compare_interpolation_methods(t_original, signal_db, method='auto'):
    """
    Универсальная функция для выбора метода интерполяции
    
    Parameters:
    -----------
    t_original : array
        Исходная временная ось
    signal_db : array
        Сигнал в dB
    method : str
        Метод интерполяции или 'auto' для автоматического выбора
    
    Returns:
    --------
    tuple : (t_interp, signal_interp_db, method_used)
    """
    
    if method == 'auto':
        # Простой автоматический выбор на основе тестов
        # В реальном коде можно добавить более сложную логику
        method = 'fft_hann'
    
    method_functions = {
        'sinc_wong': sinc_interpolation_wong,
        'fft_basic': fourier_interpolation_basic,
        'fft_hann': fourier_interpolation_hann,
        'cubic_spline': cubic_spline_interpolation,
        'linear': linear_interpolation,
        'lagrange': lambda t, s: lagrange_interpolation(t, s, 4, 3),
        'fir_filter': fir_filter_interpolation
    }
    
    if method not in method_functions:
        raise ValueError(f"Неизвестный метод: {method}")
    
    t_interp, signal_interp_db = method_functions[method](t_original, signal_db, 4)
    
    return t_interp, signal_interp_db, method

# Запуск тестирования
if __name__ == "__main__":
    print("ТЕСТИРОВАНИЕ МЕТОДОВ ИНТЕРПОЛЯЦИИ")
    print("Сравнение среднеквадратичного отклонения от теоретического sinc-сигнала")
    print("="*70)
    
    # Основной тест на идеальном sinc
    results = test_interpolation_methods()
    
    # Дополнительный тест с шумом
    noisy_results = test_with_noisy_signal()
    
    print("\n" + "="*70)
    print("ВЫВОДЫ:")
    print("="*70)
    print("1. FFT методы обычно показывают лучшую устойчивость к артефактам")
    print("2. Sinc (Wong) дает максимальную точность на идеальных сигналах") 
    print("3. FIR фильтр хорошо работает с зашумленными сигналами")
    print("4. Для большинства практических задач рекомендуется FFT + Hann")
    print("5. Линейная интерполяция - самая быстрая, но наименее точная")        
        # Идеальный sinc
        t_ideal = np.linspace(-3, 3, 30)
        ideal_sinc = np.sinc(t_ideal)
        ideal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
        
        # Зашумленный sinc
        np.random.seed(42)
        noise = np.random.normal(0, 0.05, len(ideal_sinc))
        noisy_sinc = ideal_sinc + noise
        noisy_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
        
        # Два близких пика (для проверки разрешения)
        t_peaks = np.linspace(-2, 2, 40)
        two_peaks = np.sinc(t_peaks * 2) + 0.7 * np.sinc((t_peaks - 0.5) * 2)
        peaks_db = 20 * np.log10(np.abs(two_peaks) + 1e-12)
        
        return {
            'ideal_sinc': (t_ideal, ideal_db),
            'noisy_sinc': (t_ideal, noisy_db),
            'two_peaks': (t_peaks, peaks_db)
        }
    
    def sinc_interpolation_wong(self, t_original, signal_db, interp_factor=4):
        """Sinc-интерполяция по алгоритму Wong 2005"""
        
        kernel_size = 8
        subsamples = 16
        
        def generate_coefficients(kernel_size, subsamples):
            coefficients = np.zeros((subsamples, kernel_size))
            window = np.kaiser(kernel_size, beta=2.5)

            for phase in range(subsamples):
                shift = phase / subsamples
                for k in range(kernel_size):
                    pos = (k - (kernel_size-1)//2) - shift
                    if abs(pos) < 1e-12:
                        sinc_val = 1.0
                    else:
                        sinc_val = np.sin(np.pi * pos) / (np.pi * pos)
                    coefficients[phase, k] = sinc_val * window[k]

                sum_coeff = np.sum(coefficients[phase])
                if abs(sum_coeff) > 1e-12:
                    coefficients[phase] /= sum_coeff
            return coefficients

        # Основная функция интерполяции
        coefficients_table = generate_coefficients(kernel_size, subsamples)
        signal_linear = 10 ** (signal_db / 20)
        dt = t_original[1] - t_original[0]
        
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        signal_interp = np.zeros_like(t_interp)
        half_kernel = kernel_size // 2

        for i, t_point in enumerate(t_interp):
            idx_center = int(np.floor((t_point - t_original[0]) / dt))
            fractional = (t_point - t_original[idx_center]) / dt
            
            phase_idx = int(fractional * subsamples)
            phase_idx = max(0, min(phase_idx, subsamples - 1))
            coeffs = coefficients_table[phase_idx]

            start_idx = idx_center - half_kernel
            end_idx = start_idx + kernel_size

            if start_idx < 0 or end_idx > len(signal_linear):
                signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            else:
                points = signal_linear[start_idx:end_idx]
                signal_interp[i] = np.sum(points * coeffs)

        signal_interp_db = 20 * np.log10(np.abs(signal_interp) + 1e-12)
        return t_interp, signal_interp_db
    
    def fourier_interpolation_basic(self, t_original, signal_db, interp_factor=4):
        """Базовая Фурье-интерполяция с дополнением нулями"""
        
        signal_linear = 10 ** (signal_db / 20)
        N_original = len(signal_linear)
        
        spectrum = np.fft.fft(signal_linear)
        N_new = N_original * interp_factor
        new_spectrum = np.zeros(N_new, dtype=complex)
        
        half_original = N_original // 2
        new_spectrum[:half_original] = spectrum[:half_original]
        new_spectrum[-half_original:] = spectrum[-half_original:]
        new_spectrum *= interp_factor
        
        signal_interp_linear = np.fft.ifft(new_spectrum).real
        t_interp = np.linspace(t_original[0], t_original[-1], N_new)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def fourier_interpolation_hann(self, t_original, signal_db, interp_factor=4):
        """Фурье-интерполяция с оконной функцией Hann"""
        
        signal_linear = 10 ** (signal_db / 20)
        N_original = len(signal_linear)
        
        # Оконная функция
        window = np.hanning(N_original)
        signal_windowed = signal_linear * window
        
        spectrum = np.fft.fft(signal_windowed)
        N_new = N_original * interp_factor
        new_spectrum = np.zeros(N_new, dtype=complex)
        
        half_original = N_original // 2
        new_spectrum[:half_original] = spectrum[:half_original]
        new_spectrum[-half_original:] = spectrum[-half_original:]
        new_spectrum *= interp_factor
        
        signal_interp_linear = np.fft.ifft(new_spectrum).real
        
        # Компенсация оконной функции
        window_interp = np.interp(
            np.linspace(0, N_original-1, N_new),
            np.arange(N_original),
            window
        )
        
        signal_interp_linear = np.divide(
            signal_interp_linear, 
            window_interp + 1e-12,
            out=signal_interp_linear,
            where=window_interp > 1e-6
        )
        
        t_interp = np.linspace(t_original[0], t_original[-1], N_new)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def cubic_spline_interpolation(self, t_original, signal_db, interp_factor=4):
        """Кубическая сплайн-интерполяция"""
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        
        cs = CubicSpline(t_original, signal_linear)
        signal_interp_linear = cs(t_interp)
        
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        return t_interp, signal_interp_db
    
    def linear_interpolation(self, t_original, signal_db, interp_factor=4):
        """Линейная интерполяция"""
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        
        signal_interp_linear = np.interp(t_interp, t_original, signal_linear)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def lagrange_interpolation(self, t_original, signal_db, interp_factor=4, order=3):
        """Полиномиальная интерполяция методом Лагранжа"""
        
        def lagrange_polynomial(x, x_points, y_points):
            total = 0
            n = len(x_points)
            for i in range(n):
                xi, yi = x_points[i], y_points[i]
                product = 1
                for j in range(n):
                    if i != j:
                        xj = x_points[j]
                        product *= (x - xj) / (xi - xj)
                total += yi * product
            return total
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        signal_interp_linear = np.zeros_like(t_interp)
        
        # Применяем скользящее окно Лагранжа
        for i, t_point in enumerate(t_interp):
            # Находим ближайшие точки
            distances = np.abs(t_original - t_point)
            nearest_indices = np.argsort(distances)[:order+1]
            nearest_t = t_original[nearest_indices]
            nearest_y = signal_linear[nearest_indices]
            
            signal_interp_linear[i] = lagrange_polynomial(t_point, nearest_t, nearest_y)
        
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        return t_interp, signal_interp_db
    
    def fir_filter_interpolation(self, t_original, signal_db, interp_factor=4):
        """Интерполяция через FIR фильтр (метод увеличения частоты дискретизации)"""
        
        signal_linear = 10 ** (signal_db / 20)
        
        # Создаем FIR фильтр для интерполяции
        numtaps = 31
        cutoff = 0.8 / interp_factor
        fir_coeff = firwin(numtaps, cutoff)
        
        # Увеличиваем частоту дискретизации
        upsampled = np.zeros(len(signal_linear) * interp_factor)
        upsampled[::interp_factor] = signal_linear
        
        # Применяем FIR фильтр
        signal_interp_linear = lfilter(fir_coeff, 1.0, upsampled)
        
        # Обрезаем переходный процесс
        signal_interp_linear = signal_interp_linear[numtaps//2:]
        
        t_interp = np.linspace(t_original[0], t_original[-1], len(signal_interp_linear))
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def calculate_metrics(self, t_original, signal_original_db, t_interp, signal_interp_db, method_name):
        """Вычисление метрик качества интерполяции"""
        
        # Интерполируем исходные точки для сравнения
        signal_linear_original = 10 ** (signal_original_db / 20)
        signal_linear_interp = 10 ** (signal_interp_db / 20)
        
        # Находим соответствующие точки в интерполированном сигнале
        original_points_interp = np.interp(t_original, t_interp, signal_linear_interp)
        
        # Метрики
        mse = np.mean((signal_linear_original - original_points_interp) ** 2)
        
        # Плавность (средняя производная)
        derivative = np.diff(signal_interp_db)
        smoothness = np.mean(np.abs(derivative))
        
        # Сохранение энергии
        energy_original = np.sum(signal_linear_original ** 2)
        energy_interp = np.sum(signal_linear_interp ** 2) / len(signal_linear_interp) * len(signal_linear_original)
        energy_preservation = abs(energy_original - energy_interp) / energy_original
        
        # Максимальная ошибка
        max_error = np.max(np.abs(signal_linear_original - original_points_interp))
        
        return {
            'mse': mse,
            'smoothness': smoothness,
            'energy_preservation': energy_preservation,
            'max_error': max_error
        }
    
    def test_all_methods(self, t_original, signal_db, signal_name="Сигнал"):
        """Тестирование всех методов на заданном сигнале"""
        
        interp_factor = 4
        metrics = {}
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        axes = axes.flatten()
        
        methods = [
            ('sinc_wong', self.sinc_interpolation_wong),
            ('fft_basic', self.fourier_interpolation_basic),
            ('fft_hann', self.fourier_interpolation_hann),
            ('cubic_spline', self.cubic_spline_interpolation),
            ('linear', self.linear_interpolation),
            ('lagrange', lambda t, s: self.lagrange_interpolation(t, s, interp_factor, 3)),
            ('fir_filter', self.fir_filter_interpolation)
        ]
        
        for idx, (method_key, method_func) in enumerate(methods):
            if idx >= len(axes):
                break
                
            try:
                t_interp, signal_interp = method_func(t_original, signal_db, interp_factor)
                
                # Вычисляем метрики
                metrics[method_key] = self.calculate_metrics(
                    t_original, signal_db, t_interp, signal_interp, method_key
                )
                
                # Визуализация
                ax = axes[idx]
                ax.plot(t_original, signal_db, 'bo-', label='Исходный', markersize=4, linewidth=1, alpha=0.7)
                ax.plot(t_interp, signal_interp, 'r-', label='Интерполированный', linewidth=1.5, alpha=0.8)
                ax.set_title(f'{self.methods[method_key]}', fontsize=12)
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                
                # Добавляем метрики на график
                metric_text = (f'MSE: {metrics[method_key]["mse"]:.2e}\n'
                             f'Плавность: {metrics[method_key]["smoothness"]:.3f}\n'
                             f'Энергия: {metrics[method_key]["energy_preservation"]:.2e}')
                ax.text(0.02, 0.98, metric_text, transform=ax.transAxes, fontsize=8,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
                
            except Exception as e:
                print(f"Ошибка в методе {method_key}: {e}")
                axes[idx].text(0.5, 0.5, f'Ошибка:\n{str(e)}', 
                              transform=axes[idx].transAxes, ha='center', va='center')
                axes[idx].set_title(f'{self.methods[method_key]} - ОШИБКА')
        
        # Скрываем пустые subplots
        for idx in range(len(methods), len(axes)):
            axes[idx].set_visible(False)
        
        plt.suptitle(f'Сравнение методов интерполяции: {signal_name}', fontsize=16)
        plt.tight_layout()
        plt.show()
        
        return metrics
    
    def comprehensive_test(self):
        """Комплексное тестирование на всех типах сигналов"""
        
        test_signals = self.generate_test_signals()
        all_metrics = {}
        
        for signal_name, (t, signal_db) in test_signals.items():
            print(f"\n{'='*50}")
            print(f"Тестирование: {signal_name}")
            print(f"{'='*50}")
            
            metrics = self.test_all_methods(t, signal_db, signal_name)
            all_metrics[signal_name] = metrics
            
            # Сводная таблица метрик
            self.print_metrics_table(metrics, signal_name)
        
        return all_metrics
    
    def print_metrics_table(self, metrics, signal_name):
        """Печать сводной таблицы метрик"""
        
        print(f"\nСводные метрики для {signal_name}:")
        print("-" * 80)
        print(f"{'Метод':<20} {'MSE':<12} {'Плавность':<12} {'Энергия':<12} {'Макс. ошибка':<12}")
        print("-" * 80)
        
        for method_key, metric in metrics.items():
            print(f"{self.methods[method_key]:<20} {metric['mse']:<12.2e} {metric['smoothness']:<12.3f} "
                  f"{metric['energy_preservation']:<12.2e} {metric['max_error']:<12.2e}")
    
    def find_best_method(self, metrics_dict, criterion='mse'):
        """Нахождение лучшего метода по заданному критерию"""
        
        best_methods = {}
        
        for signal_name, metrics in metrics_dict.items():
            best_score = float('inf')
            best_method = None
            
            for method_key, metric in metrics.items():
                score = metric[criterion]
                if score < best_score:
                    best_score = score
                    best_method = method_key
            
            best_methods[signal_name] = (best_method, best_score)
        
        print(f"\nЛучшие методы по критерию '{criterion}':")
        print("-" * 50)
        for signal_name, (method, score) in best_methods.items():
            print(f"{signal_name:<15} -> {self.methods[method]:<20} (значение: {score:.2e})")
        
        return best_methods

def universal_interpolation_analysis(t_original, signal_db, time_range, method='auto', interp_factor=4):
    """
    Универсальная функция анализа с автоматическим выбором или указанием метода интерполяции
    
    Parameters:
    -----------
    t_original : array
        Исходная временная ось
    signal_db : array
        Сигнал в dB
    time_range : float
        Диапазон времени (для совместимости)
    method : str
        Метод интерполяции: 'auto', 'sinc_wong', 'fft_basic', 'fft_hann', 
                           'cubic_spline', 'linear', 'lagrange', 'fir_filter'
    interp_factor : int
        Коэффициент интерполяции
    
    Returns:
    --------
    dict : Результаты анализа
    """
    
    comparator = InterpolationComparator()
    
    # Автоматический выбор метода
    if method == 'auto':
        # Простой тест для автоматического выбора
        signal_linear = 10 ** (signal_db / 20)
        noise_level = np.std(signal_linear) / np.mean(signal_linear)
        
        if noise_level > 0.1:  # Зашумленный сигнал
            method = 'fft_hann'
        elif len(t_original) < 20:  # Мало точек
            method = 'cubic_spline'
        else:  # Стандартный случай
            method = 'sinc_wong'
    
    # Применение выбранного метода
    method_functions = {
        'sinc_wong': comparator.sinc_interpolation_wong,
        'fft_basic': comparator.fourier_interpolation_basic,
        'fft_hann': comparator.fourier_interpolation_hann,
        'cubic_spline': comparator.cubic_spline_interpolation,
        'linear': comparator.linear_interpolation,
        'lagrange': lambda t, s: comparator.lagrange_interpolation(t, s, interp_factor, 3),
        'fir_filter': comparator.fir_filter_interpolation
    }
    
    if method not in method_functions:
        raise ValueError(f"Неизвестный метод: {method}. Доступные: {list(method_functions.keys())}")
    
    t_interp, sinc_interp = method_functions[method](t_original, signal_db, interp_factor)
    
    # Анализ характеристик (упрощенный)
    def find_main_lobe_width(t, signal_db):
        center_idx = np.argmax(signal_db)
        threshold = signal_db[center_idx] - 3
        
        # Поиск точек пересечения -3 dB
        left_idx = center_idx
        while left_idx > 0 and signal_db[left_idx] > threshold:
            left_idx -= 1
        
        right_idx = center_idx
        while right_idx < len(signal_db) - 1 and signal_db[right_idx] > threshold:
            right_idx += 1
        
        if left_idx > 0 and right_idx < len(signal_db) - 1:
            # Линейная интерполяция
            wl = t[left_idx] + (t[left_idx+1] - t[left_idx]) * (threshold - signal_db[left_idx]) / (signal_db[left_idx+1] - signal_db[left_idx])
            wr = t[right_idx-1] + (t[right_idx] - t[right_idx-1]) * (threshold - signal_db[right_idx-1]) / (signal_db[right_idx] - signal_db[right_idx-1])
            width = wr - wl
            return wl, wr, width
        return None, None, 0
    
    wl, wr, width = find_main_lobe_width(t_interp, sinc_interp)
    
    # Упрощенный расчет УБЛ
    def calculate_sidelobe_levels(t, signal_db):
        signal_linear = 10 ** (signal_db / 20)
        main_lobe_max = np.max(signal_linear)
        main_lobe_idx = np.argmax(signal_linear)
        
        # Боковые лепестки - все кроме области вокруг главного лепестка
        sidelobe_mask = np.ones(len(signal_linear), dtype=bool)
        sidelobe_mask[max(0, main_lobe_idx-5):min(len(signal_linear), main_lobe_idx+6)] = False
        
        if np.any(sidelobe_mask):
            max_sidelobe = np.max(signal_linear[sidelobe_mask])
            classical_pslr = 20 * np.log10(max_sidelobe / main_lobe_max)
            
            # Упрощенный интегральный УБЛ
            power_main = np.sum(signal_linear[~sidelobe_mask] ** 2)
            power_sidelobes = np.sum(signal_linear[sidelobe_mask] ** 2)
            integral_pslr = 10 * np.log10(power_sidelobes / power_main)
        else:
            classical_pslr = integral_pslr = -80
        
        return classical_pslr, integral_pslr
    
    classica
