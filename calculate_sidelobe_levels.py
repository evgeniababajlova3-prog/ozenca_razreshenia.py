import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import find_peaks



def calculate_sidelobe_levels(t, signal_db):

    # Конвертируем в линейную область для анализа
    signal_linear = 10 ** (signal_db / 20)

    # Находим главный максимум
    main_peak_idx = np.argmax(signal_linear)
    main_peak_val = signal_linear[main_peak_idx]

    # Находим все локальные минимумы (нули sinc функции)
    minima_indices, _ = find_peaks(-signal_linear, height=-0.1, distance=5)

    if len(minima_indices) < 2:
        print("Не удалось найти достаточное количество нулей функции")
        return -80, -80

    # Разделяем минимумы на левые и правые относительно главного максимума
    left_minima = minima_indices[minima_indices < main_peak_idx]
    right_minima = minima_indices[minima_indices > main_peak_idx]

    if len(left_minima) == 0 or len(right_minima) == 0:
        print("Не удалось найти нули по обе стороны от главного максимума")
        return -80, -80

    # Берем ближайшие нули слева и справа от главного максимума
    left_zero_idx = left_minima[-1]  # последний нуль слева от максимума
    right_zero_idx = right_minima[0]  # первый нуль справа от максимума

    # Главный лепесток - между этими нулями
    main_lobe_indices = range(left_zero_idx, right_zero_idx + 1)
    main_lobe = signal_linear[main_lobe_indices]
    t_main = t[main_lobe_indices]

    # Боковые лепестки - все остальное
    sidelobe_indices = np.concatenate([
        np.arange(0, left_zero_idx),
        np.arange(right_zero_idx + 1, len(signal_linear))
    ])
    sidelobes = signal_linear[sidelobe_indices]
    t_sidelobes = t[sidelobe_indices]

    # Классический УБЛ - отношение максимального бокового лепестка к главному
    if len(sidelobes) > 0:
        max_sidelobe = np.max(sidelobes)
        classical_pslr_db = 20 * np.log10(max_sidelobe / main_peak_val)
    else:
        classical_pslr_db = -80

    # Интегральный УБЛ - отношение мощностей
    if len(main_lobe) > 0 and len(sidelobes) > 0:
        power_main = trapezoid(main_lobe ** 2, t_main)
        power_sidelobes = trapezoid(sidelobes ** 2, t_sidelobes)
        integral_pslr_db = 10 * np.log10(power_sidelobes / power_main)
    else:
        integral_pslr_db = -80

    return classical_pslr_db, integral_pslr_db