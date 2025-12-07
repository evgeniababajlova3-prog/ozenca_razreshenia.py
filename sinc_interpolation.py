from generate_coefficients import generate_coefficients
import numpy as np

def sinc_interpolation(t_original, signal_db, kernel_size, subsamples, f_discr, interp_factor):

    # Генерируем коэффициенты
    coefficients_table = generate_coefficients(kernel_size, subsamples)

    # Создаем интерполированную временную сетку
    t_interp = np.linspace(t_original[0], t_original[-1], f_discr * interp_factor)
    signal_interp = np.zeros_like(t_interp)

    half_kernel = kernel_size // 2

    # Преобразуем сигнал в линейную область для интерполяции
    signal_linear = 10 ** (signal_db / 20)

    # Шаг дискретизации исходного сигнала
    dt = t_original[1] - t_original[0]

    for i, t_point in enumerate(t_interp):
    # Находим ближайший индекс слева в исходной сетке
        idx_center = int(np.floor((t_point - t_original[0]) / dt))

    # Проверяем границы
        if idx_center < half_kernel or idx_center >= len(t_original) - half_kernel:
        # Для граничных точек используем линейную интерполяцию
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            continue

    # Определяем дробную часть сдвига
        fractional = (t_point - t_original[idx_center]) / dt

    # Индекс в таблице коэффициентов
        phase_idx = int(fractional * subsamples)
        phase_idx = max(0, min(phase_idx, subsamples - 1))

    # Берем коэффициенты для данной фазы
        coeffs = coefficients_table[phase_idx]

    # Индексы исходных отсчетов для свертки
        start_idx = idx_center - half_kernel
        end_idx = start_idx + kernel_size

    # Проверяем границы еще раз
        if start_idx < 0 or end_idx > len(signal_linear):
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            continue

    # Выполняем свертку
        points = signal_linear[start_idx:end_idx]
        signal_interp[i] = np.sum(points * coeffs)

    # Преобразуем обратно
    signal_interp_db = 20 * np.log10(np.abs(signal_interp))

    return t_interp, signal_interp_db