import matplotlib.pyplot as plt
import numpy as np

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
        ax2.axhline(y=-3, color='r', linestyle='--', linewidth=2, label=f'Уровень 3 дБ')


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
        ax2.plot(wl, -3, 'mo', markersize=8)

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
        ax2.plot(wr, -3, 'mo', markersize=8, label='Точки пересечения')

    # Показываем ширину
        ax2.plot([wl, wr], [-3, -3], 'm-', linewidth=3, alpha=0.7, label=f'Ширина: {width:.4f}')

        ax2.set_xlabel('t')
        ax2.set_ylabel('Амплитуда (дБ)')
        ax2.set_title('Детальный вид точек пересечения')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(-3 - 1, -3 + 1)
    else:
        ax2.text(0.5, 0.5, 'Не удалось определить\nграницы лепестка',
             horizontalalignment='center', verticalalignment='center',
             transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Детальный вид точек пересечения')

    plt.tight_layout()
    plt.show()