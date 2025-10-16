import numpy as np
from scipy.integrate import trapezoid

def calculate_sidelobe_levels(t, signal_db, sinc_width):

    # Конвертируем в линейную область для анализа мощности
    signal_linear = 10 ** (signal_db / 20)

    # Маска для главного лепестка
    main_lobe_mask = (t >= -sinc_width) & (t <= sinc_width)

    # Главный лепесток
    main_lobe = signal_linear[main_lobe_mask]
    t_main = t[main_lobe_mask]

    # Боковые лепестки (все, кроме главного)

    sidelobes = signal_linear[~main_lobe_mask]
    t_sidelobes = t[~main_lobe_mask]

    if len(main_lobe) == 0 or len(sidelobes) == 0:
        return -80, -80

    #Классический УБЛ - отношение максимального бокового лепестка к главному
    max_main = np.max(main_lobe)
    max_sidelobe = np.max(sidelobes)
    classical_pslr_db = 20 * np.log10(max_sidelobe / max_main) if max_sidelobe > 0 else -80

    #Интегральный УБЛ - отношение мощностей
    power_main = trapezoid(main_lobe ** 2, t_main)
    power_sidelobes = trapezoid(sidelobes ** 2, t_sidelobes)
    integral_pslr_db = 10 * np.log10(power_sidelobes / power_main) if power_sidelobes > 0 else -80

    return classical_pslr_db, integral_pslr_db