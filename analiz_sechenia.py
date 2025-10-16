from sinc_interpolation import sinc_interpolation
from theoretical_sinc_width import theoretical_sinc_width
from find_main_lobe_width import find_main_lobe_width
from calculate_sidelobe_levels import calculate_sidelobe_levels
from plot_results import plot_results

# Параметры сигнала
F_DISCR = 200  # Количество точек дискретизации [рекомендуется: 200-2000]
SINC_WIDTH = 1.0  # Ширина sinc-функции [рекомендуется: 0.5-2.0]
TIME_RANGE = 10.0  # Диапазон времени [-TIME_RANGE, TIME_RANGE]
# Параметры интерполяции
KERNEL_SIZE = 8  # Размер ядра [рекомендуется: 6-12, должно быть четным]
SUBSAMPLES = 16  # Количество подвыборок [рекомендуется: 8-32]
INTERP_FACTOR = 10  # Коэффициент интерполяции [рекомендуется: 2-10]
# Параметры анализа
THRESHOLD_DB = -3  # Уровень для определения ширины лепестка [дБ]
MAIN_LOBE_SEARCH_RADIUS = 2.0  # Радиус поиска главного лепестка
THEORETICAL_CLASSICAL_PSLR = -13.26  # Классический УБЛ [дБ]
THEORETICAL_INTEGRAL_PSLR = -10.96   # Интегральный УБЛ [дБ]

#Основная функция
def analiz_sechenia(t_original, sinc_db, sinc_linear):

    # Исходные данные
    print("\nИСХОДНЫЕ ДАННЫЕ:")
    print(f"Количество точек: {F_DISCR}")
    print(f"Ширина sinc-функции: {SINC_WIDTH}")
    print(f"Размер ядра: {KERNEL_SIZE}")
    print(f"Количество фаз: {SUBSAMPLES}")
    print(f"Коэффициент интерполяции: {INTERP_FACTOR}")

    #Sinc-интерполяция
    t_interp, sinc_interp = sinc_interpolation(t_original, sinc_db, KERNEL_SIZE, SUBSAMPLES, TIME_RANGE, F_DISCR, INTERP_FACTOR)

    #Нахождение ширины главного лепестка
    wl, wr, width, left_points, right_points = find_main_lobe_width(t_interp, sinc_interp, THRESHOLD_DB)

    # Теоретическое значение
    theoretical_width = theoretical_sinc_width(SINC_WIDTH)

    print(f"\nРЕЗУЛЬТАТЫ ИЗМЕРЕНИЙ:")
    print(f"\nТеоретическое разрешение: {theoretical_width:.6f}")
    if width is not None:
        print(f"Измеренное разрешение:  {width:.6f}")
        error_width = abs(width - theoretical_width) / theoretical_width * 100
        print(f"Погрешность измерения разрешения:  {error_width:.2f}%")
    else:
        print("Не удалось определить ширину лепестка")
        width = 0

    if wl is not None and wr is not None:
        classical_pslr, integral_pslr = calculate_sidelobe_levels(t_interp, sinc_interp, SINC_WIDTH)
        print(f"\nТеоретический классический УБЛ: {THEORETICAL_CLASSICAL_PSLR:.2f} дБ")
        print(f"Измеренный классический УБЛ:   {classical_pslr:.2f} дБ")
        error_classical_pslr = abs(classical_pslr - THEORETICAL_CLASSICAL_PSLR)
        print(f"Погрешность измерения классического УБЛ:  {error_classical_pslr:.2f}дБ")
        print(f"\nТеоретический интегральный УБЛ: {THEORETICAL_INTEGRAL_PSLR:.2f} дБ")
        print(f"Измеренный интегральный УБЛ:   {integral_pslr:.2f} дБ")
        error_integral_pslr = abs(integral_pslr - THEORETICAL_INTEGRAL_PSLR)
        print(f"Погрешность измерения интегрального УБЛ:  {error_integral_pslr:.2f}дБ")
    else:
        classical_pslr = integral_pslr = -80
        print("Не удалось вычислить УБЛ (не определены границы главного лепестка)")

    # 5. Визуализация
    plot_results(t_original, sinc_db, t_interp, sinc_interp, wl, wr, width, left_points, right_points, THRESHOLD_DB)


    return {
        'theoretical_width': theoretical_width,
        'measured_width': width,
        'classical_pslr': classical_pslr,
        'integral_pslr': integral_pslr
    }

