import numpy as np

def find_main_lobe_width(t, signal_db, threshold_db):
    center_idx = np.argmax(signal_db)

    # Ищем точки слева от центра
    left_points = []
    for i in range(center_idx, 0, -1):
        if (signal_db[i - 1] - threshold_db) * (signal_db[i] - threshold_db) <= 0:
            # Нашли интервал пересечения, берем две точки для интерполяции
            left_points = [(t[i - 1], signal_db[i - 1]), (t[i], signal_db[i])]
            break

    # Ищем точки справа от центра
    right_points = []
    for i in range(center_idx, len(signal_db) - 1):
        if (signal_db[i] - threshold_db) * (signal_db[i + 1] - threshold_db) <= 0:
            # Нашли интервал пересечения, берем две точки для интерполяции
            right_points = [(t[i], signal_db[i]), (t[i + 1], signal_db[i + 1])]
            break

    if not left_points or not right_points:
        return None, None, None, None, None

    # Вычисляем точки пересечения через линейную интерполяцию
    # Для левой стороны
    x1, y1 = left_points[0]
    x2, y2 = left_points[1]
    wl = x1 + (x2 - x1) * (threshold_db - y1) / (y2 - y1)

    # Для правой стороны
    x1, y1 = right_points[0]
    x2, y2 = right_points[1]
    wr = x1 + (x2 - x1) * (threshold_db - y1) / (y2 - y1)

    width = wr - wl

    return wl, wr, width, left_points, right_points