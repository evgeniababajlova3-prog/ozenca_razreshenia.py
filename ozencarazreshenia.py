import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

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

#Генерация нормированного sinc-сигнала с заданной шириной
def generate_sinc_signal():
    t = np.linspace(-TIME_RANGE, TIME_RANGE, F_DISCR)
    sinc_linear = np.sinc(t / SINC_WIDTH)
    sinc_norm = np.abs(sinc_linear) / np.max(np.abs(sinc_linear))
    sinc_db = 20 * np.log10(sinc_norm + 1e-12)
    sinc_db = np.clip(sinc_db, -80, 0)

    return t, sinc_db, sinc_norm

#Генерация коэффициентов
def generate_coefficients(kernel_size, subsamples):
    coefficients = np.zeros((subsamples, kernel_size))

    # Создаем окно Кайзера для уменьшения эффекта Гиббса
    window = np.kaiser(kernel_size, beta=2.5)

    for phase in range(subsamples):
        # Фаза интерполяции (0 = целое число, subsamples-1 = почти следующее целое)
        shift = phase / subsamples

        for k in range(kernel_size):
            # Позиция точки относительно центра интерполяции
            # Центр ядра находится между точками при kernel_size четном
            pos = (k - kernel_size // 2 + 0.5) - shift

            # Вычисляем sinc функцию
            if abs(pos) < 1e-12:
                sinc_val = 1.0
            else:
                sinc_val = np.sin(np.pi * pos) / (np.pi * pos)

            # Применяем оконную функцию
            coefficients[phase, k] = sinc_val * window[k]

    # Нормализация коэффициентов (сумма = 1)
    for phase in range(subsamples):
        sum_coeff = np.sum(coefficients[phase])
        if abs(sum_coeff) > 1e-12:
            coefficients[phase] /= sum_coeff

    return coefficients

#Sinc-интерполяция
def wong_sinc_interpolation(t_original, signal_db):

    # Генерируем коэффициенты
    coefficients_table = generate_coefficients(KERNEL_SIZE, SUBSAMPLES)

    # Создаем интерполированную временную сетку
    t_interp = np.linspace(-TIME_RANGE, TIME_RANGE, F_DISCR * INTERP_FACTOR)
    signal_interp = np.zeros_like(t_interp)

    half_kernel = KERNEL_SIZE // 2

    # Преобразуем сигнал в линейную область для интерполяции
    signal_linear = 10 ** (signal_db / 20)

    # Шаг дискретизации исходного сигнала
    dt = t_original[1] - t_original[0]

    for i, t_point in enumerate(t_interp):
    # Находим ближайший индекс слева в исходной сетке
        idx_center = int((t_point - t_original[0]) / dt)

    # Проверяем границы
        if idx_center < half_kernel or idx_center >= len(t_original) - half_kernel:
        # Для граничных точек используем линейную интерполяцию
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            continue

    # Определяем дробную часть сдвига
        fractional = (t_point - t_original[idx_center]) / dt

    # Индекс в таблице коэффициентов
        phase_idx = int(fractional * SUBSAMPLES)
        phase_idx = max(0, min(phase_idx, SUBSAMPLES - 1))

    # Берем коэффициенты для данной фазы
        coeffs = coefficients_table[phase_idx]

    # Индексы исходных отсчетов для свертки
        start_idx = idx_center - half_kernel
        end_idx = start_idx + KERNEL_SIZE

    # Проверяем границы еще раз
        if start_idx < 0 or end_idx > len(signal_linear):
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            continue

    # Выполняем свертку
        points = signal_linear[start_idx:end_idx]
        signal_interp[i] = np.sum(points * coeffs)

    # Преобразуем обратно в дБ
    signal_interp_db = 20 * np.log10(np.abs(signal_interp) + 1e-12)
    signal_interp_db = np.clip(signal_interp_db, -80, 0)

    return t_interp, signal_interp_db


#Теоретическая ширина главного лепестка sinc-функции на уровне -3 дБ
def theoretical_sinc_width(sinc_width):
    return 0.8859 * sinc_width

#Нахождение ширины главного лепестка на уровне -3 дБ
def find_main_lobe_width(t, signal_db):
    center_idx = np.argmax(signal_db)

    # Ищем точки слева от центра
    left_points = []
    for i in range(center_idx, 0, -1):
        if (signal_db[i - 1] - THRESHOLD_DB) * (signal_db[i] - THRESHOLD_DB) <= 0:
            # Нашли интервал пересечения, берем две точки для интерполяции
            left_points = [(t[i - 1], signal_db[i - 1]), (t[i], signal_db[i])]
            break

    # Ищем точки справа от центра
    right_points = []
    for i in range(center_idx, len(signal_db) - 1):
        if (signal_db[i] - THRESHOLD_DB) * (signal_db[i + 1] - THRESHOLD_DB) <= 0:
            # Нашли интервал пересечения, берем две точки для интерполяции
            right_points = [(t[i], signal_db[i]), (t[i + 1], signal_db[i + 1])]
            break

    if not left_points or not right_points:
        return None, None, None, None, None

    # Вычисляем точки пересечения через линейную интерполяцию
    # Для левой стороны
    x1, y1 = left_points[0]
    x2, y2 = left_points[1]
    wl = x1 + (x2 - x1) * (THRESHOLD_DB - y1) / (y2 - y1)

    # Для правой стороны
    x1, y1 = right_points[0]
    x2, y2 = right_points[1]
    wr = x1 + (x2 - x1) * (THRESHOLD_DB - y1) / (y2 - y1)

    width = wr - wl

    return wl, wr, width, left_points, right_points

#Вычисление УБЛ
def calculate_sidelobe_levels(t, signal_db, SINC_WIDTH):

    # Конвертируем в линейную область для анализа мощности
    signal_linear = 10 ** (signal_db / 20)

    # Маска для главного лепестка
    main_lobe_mask = (t >= -SINC_WIDTH) & (t <= SINC_WIDTH)

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

#Визуализация
def plot_results(t_original, sinc_db, t_interp, sinc_interp, wl, wr, width, left_points, right_points):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # 1. Исходный и интерполированный сигнал
    ax1.plot(t_original, sinc_db, 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
    ax1.plot(t_interp, sinc_interp, 'g-', linewidth=1, label='Интерполированный сигнал', alpha=0.7)
    ax1.set_xlabel('t')
    ax1.set_ylabel('Амплитуда (дБ)')
    ax1.set_title('Сравнение исходного и интерполированного сигнала')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-50, 5)
    ax1.set_xlim(-10, 10)

    # 2. Детальный вид точек пересечения
    if wl is not None and wr is not None and left_points and right_points:
        zoom_margin = width * 0.5
        mask = (t_interp >= wl - zoom_margin) & (t_interp <= wr + zoom_margin)
        t_zoom = t_interp[mask]
        sinc_zoom = sinc_interp[mask]

        # Рисуем интерполированный сигнал
        ax2.plot(t_zoom, sinc_zoom, 'g-', linewidth=2, alpha=0.7, label='Интерполированный сигнал')
        ax2.axhline(y=THRESHOLD_DB, color='r', linestyle='--', linewidth=2, label=f'Уровень {THRESHOLD_DB} дБ')


    # Левая сторона - точки и прямая
        left_t = [p[0] for p in left_points]
        left_y = [p[1] for p in left_points]

    # Вычисляем параметры прямой для левой стороны
        k_left = (left_y[1] - left_y[0]) / (left_t[1] - left_t[0])
        b_left = left_y[0] - k_left * left_t[0]

    # Продлеваем прямую за пределы точек
        t_left_extended = np.array([left_t[0] - 0.1, left_t[1] + 0.1])
        y_left_extended = k_left * t_left_extended + b_left

        ax2.plot(t_left_extended, y_left_extended, 'b--', linewidth=2)
        ax2.plot(left_t, left_y, 'bo', markersize=6)
        ax2.plot(wl, THRESHOLD_DB, 'mo', markersize=8)

    # Правая сторона - точки и прямая
        right_t = [p[0] for p in right_points]
        right_y = [p[1] for p in right_points]

    # Вычисляем параметры прямой для правой стороны
        k_right = (right_y[1] - right_y[0]) / (right_t[1] - right_t[0])
        b_right = right_y[0] - k_right * right_t[0]

    #Продлеваем прямую за пределы точек
        t_right_extended = np.array([right_t[0] - 0.1, right_t[1] + 0.1])
        y_right_extended = k_right * t_right_extended + b_right

        ax2.plot(t_right_extended, y_right_extended, 'b--', linewidth=2, label='Интерполяционные прямые')
        ax2.plot(right_t, right_y, 'bo', markersize=6, label='Точки вблизи точки пересечения')
        ax2.plot(wr, THRESHOLD_DB, 'mo', markersize=8, label='Точки пересечения')

    # Показываем ширину
        ax2.plot([wl, wr], [THRESHOLD_DB, THRESHOLD_DB], 'm-', linewidth=3, alpha=0.7, label=f'Ширина: {width:.4f}')

        ax2.set_xlabel('t')
        ax2.set_ylabel('Амплитуда (дБ)')
        ax2.set_title('Детальный вид точек пересечения')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(THRESHOLD_DB - 1, THRESHOLD_DB + 1)
    else:
        ax2.text(0.5, 0.5, 'Не удалось определить\nграницы лепестка',
             horizontalalignment='center', verticalalignment='center',
             transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Детальный вид точек пересечения')

    plt.tight_layout()
    plt.show()

#Основная функция
def main():

    # Исходные данные
    print("\nИСХОДНЫЕ ДАННЫЕ:")
    print(f"Количество точек: {F_DISCR}")
    print(f"Ширина sinc-функции: {SINC_WIDTH}")
    print(f"Размер ядра: {KERNEL_SIZE}")
    print(f"Количество фаз: {SUBSAMPLES}")
    print(f"Коэффициент интерполяции: {INTERP_FACTOR}")

    #Формирование исходных данных
    t_original, sinc_db, sinc_linear = generate_sinc_signal()

    #Sinc-интерполяция
    t_interp, sinc_interp = wong_sinc_interpolation(t_original, sinc_db)

    #Нахождение ширины главного лепестка
    wl, wr, width, left_points, right_points = find_main_lobe_width(t_interp, sinc_interp)

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
    plot_results(t_original, sinc_db, t_interp, sinc_interp, wl, wr, width, left_points, right_points)


    return {
        'theoretical_width': theoretical_width,
        'measured_width': width,
        'classical_pslr': classical_pslr,
        'integral_pslr': integral_pslr
    }

main()
