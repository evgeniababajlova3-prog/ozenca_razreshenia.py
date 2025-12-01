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