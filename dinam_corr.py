def calculate_incidence_angle_matrix(height, width, first_delay, last_delay, platform_height):
    incidence_angles = np.zeros((height, width))
    for j in range(width): # Цикл по строкам дальности
        # 1. Расчет наклонной дальности для j-го отсчета
        range_delay = first_delay + (last_delay - first_delay) * (j / width)
        slant_range = range_delay * 3e8 / 2 # м

        # 2. Расчет угла падения (в радианах)
        # Для платформы с постоянной высотой: sin(θ) = высота / наклонная дальность
        # Учитываем случай, когда slant_range меньше высоты (нефизично для дальних отсчетов)
        sin_theta = platform_height / slant_range if slant_range > platform_height else 1.0
        theta_rad = np.arcsin(sin_theta)

        incidence_angles[:, j] = np.degrees(theta_rad) # Сохраняем в градусах для всех строк азимута
    return incidence_angles

def apply_basic_radiometric_correction(complex_image, incidence_angles_deg):
    """
    Применяет поправку за наклонную дальность: σ⁰ = β⁰ / sin(θ)
    """
    # 1. Вычисляем исходную мощность (β⁰ пропорциональна |DN|^2)
    power_original = np.abs(complex_image) ** 2

    # 2. Рассчитываем поправочный коэффициент. Добавляем малый эпсилон для избежания деления на ноль.
    # Углы преобразуем обратно в радианы для вычисления синуса.
    incidence_angles_rad = np.radians(incidence_angles_deg)
    correction_factor = np.sin(incidence_angles_rad)
    correction_factor = np.where(correction_factor < 1e-12, 1e-12, correction_factor) # Защита от нуля

    # 3. Применяем коррекцию к мощности
    power_corrected = power_original / correction_factor

    # 4. Восстанавливаем комплексное изображение с новой амплитудой и сохраненной фазой
    phase_original = np.angle(complex_image)
    amplitude_corrected = np.sqrt(power_corrected)
    complex_corrected = amplitude_corrected * np.exp(1j * phase_original)

    return complex_corrected, incidence_angles_deg, correction_factor

def main():
    # ... (загрузка данных: radar_image_complex, radar_image = load_radar_image_from_hdf5(...))

    # БЛОК РАДИОМЕТРИЧЕСКОЙ КОРРЕКЦИИ (НОВЫЙ)
    print("Применение радиометрической коррекции...")
    radar_image_complex_corr, incidence_angles, corr_factor = apply_radiometric_correction(
        radar_image_complex, JSON_FILE_PATH
    )
    # Далее использовать radar_image_complex_corr вместо radar_image_complex
    corrected_amplitude_image = np.abs(radar_image_complex_corr)

    # ... (далее: обнаружение целей, анализ сечений и т.д. на основе скорректированного изображения)
    detected_peaks = find_targets(radar_image_complex_corr, MIN_DISTANCE, THRESHOLD_OFFSET_DB)

    import numpy as np

    def sinc_interpolation(t_original, signal_db, kernel_size=8, subsamples=16,
                           time_range=None, f_discr=None, interp_factor=10, beta=2.5):
        """
        Sinc-интерполяция через FFT (частотная область).
        Совместимый интерфейс с оригинальной sinc_interpolation.

        Параметры kernel_size, subsamples, beta игнорируются (оставлены для совместимости).
        Параметры time_range, f_discr также игнорируются - временная ось строится автоматически.
        """

        # 1. Преобразование из дБ в линейную шкалу
        signal_linear = 10 ** (signal_db / 20.0)

        N = len(signal_linear)
        L = interp_factor  # Коэффициент интерполяции
        N_new = N * L

        # 2. ДПФ исходного сигнала
        spectrum = np.fft.fft(signal_linear)

        # 3. Создание нового спектра с нулями
        new_spectrum = np.zeros(N_new, dtype=complex)

        # 4. Копирование спектральных компонент с учетом четности
        half_N = N // 2

        # Низкие частоты в начало
        new_spectrum[:half_N] = spectrum[:half_N]

        # Высокие частоты в конец (для сохранения симметрии)
        new_spectrum[-half_N:] = spectrum[-half_N:]

        # Если N нечетное - копируем центральный компонент
        if N % 2 != 0:
            new_spectrum[half_N] = spectrum[half_N]

        # 5. Обратное ДПФ и масштабирование
        interpolated_linear = np.real(np.fft.ifft(new_spectrum)) * L

        # 6. Преобразование обратно в дБ
        interpolated_db = 20 * np.log10(np.abs(interpolated_linear) + 1e-12)

        # 7. Создание новой временной оси (на том же интервале)
        t_interp = np.linspace(t_original[0], t_original[-1], N_new)

        return t_interp, interpolated_db