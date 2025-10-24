from sinc_interpolation import sinc_interpolation
from find_main_lobe_width import find_main_lobe_width
from calculate_sidelobe_levels import calculate_sidelobe_levels
from plot_results import plot_results

# Параметры сигнала
F_DISCR = 200  # Количество точек дискретизации [рекомендуется: 200-2000]
TIME_RANGE = 10.0  # Диапазон времени [-TIME_RANGE, TIME_RANGE]
# Параметры интерполяции
KERNEL_SIZE = 8  # Размер ядра [рекомендуется: 6-12, должно быть четным]
SUBSAMPLES = 16  # Количество подвыборок [рекомендуется: 8-32]
INTERP_FACTOR = 10  # Коэффициент интерполяции [рекомендуется: 2-10]
# Параметры анализа
THEORETICAL_CLASSICAL_PSLR = -13.26  # Классический УБЛ [дБ]
THEORETICAL_INTEGRAL_PSLR = -10.96   # Интегральный УБЛ [дБ]

#Основная функция
def analiz_sechenia(t_original, sinc_db):

    #Sinc-интерполяция
    t_interp, sinc_interp = sinc_interpolation(t_original, sinc_db, KERNEL_SIZE, SUBSAMPLES, TIME_RANGE, F_DISCR, INTERP_FACTOR)

    #Нахождение ширины главного лепестка
    wl, wr, width, left_points, right_points = find_main_lobe_width(t_interp, sinc_interp)

    print(f"\nРЕЗУЛЬТАТЫ ИЗМЕРЕНИЙ:")
    if width is not None:
        print(f"Ширина главного лепестка:  {width:.3f}")

    if wl is not None and wr is not None:
        classical_pslr, integral_pslr = calculate_sidelobe_levels(t_interp, sinc_interp)
        print(f"Максимальный УБЛ:   {classical_pslr:.2f} дБ")
        print(f"Интегральный УБЛ:   {integral_pslr:.2f} дБ")
    else:
        classical_pslr = integral_pslr = -80
        print("Не удалось вычислить УБЛ (не определены границы главного лепестка)")

    # 5. Визуализация
    plot_results(t_original, sinc_db, t_interp, sinc_interp, wl, wr, width, classical_pslr)

    return {
        'measured_width': width,
        'classical_pslr': classical_pslr,
        'integral_pslr': integral_pslr
    }

