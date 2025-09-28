import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

#Формирование исходных данных
f_discr = int(input("Кол-во точек на отсчет:"))
t = np.linspace(-10, 10, f_discr)
sinc_ish = np.sinc(t)
sinc_abs = np.abs(sinc_ish)
sinc_norm = sinc_abs / np.max(sinc_abs)
sinc_db = 20 * np.log10(sinc_norm)
sinc_db = np.nan_to_num(sinc_db, nan=-80, neginf=-80)

#Визуализация исходного сигнала
plt.figure(figsize=(12, 6))
plt.plot(t, sinc_db, label='Исходный сигнал (дБ)')
plt.axhline(y=-3, color='r', linestyle='--', alpha=0.7, label='Уровень -3 дБ')
plt.xlabel('Время')
plt.ylabel('Амплитуда (дБ)')
plt.title('Исходный нормированный сигнал sinc в дБ')
plt.legend()
plt.grid(True)
plt.ylim(-50, 5)
plt.show()

#Кубическая интерполяция
def cubic_interpolate():
    func_interp = interp1d(t, sinc_db, kind='cubic')
    k = int(input("Коэффициент кубической интерполяции:"))
    t_interp = np.linspace(-10, 10, f_discr*k)
    sinc_interp = func_interp(t_interp)
    return t_interp, sinc_interp
t_interp_1, sinc_interp_1 = cubic_interpolate()

#Визуализация кубически интерполированного сигнала
plt.figure(figsize=(12, 6))
plt.plot(t_interp_1, sinc_interp_1, label='Кубически интерполированный сигнал (дБ)')
plt.axhline(y=-3, color='r', linestyle='--', alpha=0.7, label='Уровень -3 дБ')
plt.xlabel('Время')
plt.ylabel('Амплитуда (дБ)')
plt.title('Кубически интерполированный сигнал sinc в дБ')
plt.legend()
plt.grid(True)
plt.ylim(-50, 5)
plt.show()

#Интерполяция по таблице коэффициентов
def wong_interpolate():
    coefficients_table = np.array([
        [-0.003, 0.010, -0.024, 0.062, 0.993, -0.054, 0.021, -0.009],
        [-0.007, 0.021, -0.049, 0.131, 0.973, -0.098, 0.040, -0.017],
        [-0.012, 0.032, -0.075, 0.207, 0.941, -0.134, 0.055, -0.023],
        [-0.016, 0.043, -0.101, 0.287, 0.896, -0.160, 0.066, -0.027],
        [-0.020, 0.054, -0.125, 0.371, 0.841, -0.176, 0.074, -0.030],
        [-0.024, 0.063, -0.147, 0.457, 0.776, -0.185, 0.078, -0.031],
        [-0.027, 0.071, -0.165, 0.542, 0.703, -0.185, 0.079, -0.031],
        [-0.030, 0.076, -0.178, 0.625, 0.625, -0.178, 0.076, -0.030],
        [-0.031, 0.079, -0.185, 0.703, 0.542, -0.165, 0.071, -0.027],
        [-0.031, 0.078, -0.185, 0.776, 0.457, -0.147, 0.063, -0.024],
        [-0.030, 0.074, -0.176, 0.841, 0.371, -0.125, 0.054, -0.020],
        [-0.027, 0.066, -0.160, 0.896, 0.287, -0.101, 0.043, -0.016],
        [-0.023, 0.055, -0.134, 0.941, 0.207, -0.075, 0.032, -0.012],
        [-0.017, 0.040, -0.098, 0.973, 0.131, -0.049, 0.021, -0.007],
        [-0.009, 0.021, -0.054, 0.993, 0.062, -0.024, 0.010, -0.003],
        [0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000]
    ])
    faza = 16
    k = int(input("Коэффициент интерполяции для интерполяции по таблице:"))
    t_interp = np.linspace(-10, 10, f_discr*k)
    sinc_interp = np.zeros_like(t_interp)
    for i, t_znach in enumerate(t_interp):
        idx = np.searchsorted(t, t_znach )
        if idx < 3 or idx > len(t) - 5:
            sinc_interp[i] = np.interp(t_znach , t, sinc_db)
            continue
        fraction = (t_znach  - t[idx - 1]) / (t[idx] - t[idx - 1])
        table_idx = int(fraction * (faza - 1))
        table_idx = max(0, min(table_idx, faza - 1))
        coeffs = coefficients_table[table_idx]
        points = sinc_db[idx - 3:idx + 5]
        if len(points) == 8:
            sinc_interp[i] = np.dot(points, coeffs)
        else:
            sinc_interp[i] = np.interp(t_znach, t, sinc_db)
    return t_interp, sinc_interp
t_interp_2, sinc_interp_2 = wong_interpolate()

#Визуализация сигнала, интерполированного по таблице
plt.figure(figsize=(12, 6))
plt.plot(t_interp_2, sinc_interp_2, label='Интерполированный с помощью таблицы сигнал (дБ)')
plt.axhline(y=-3, color='r', linestyle='--', alpha=0.7, label='Уровень -3 дБ')
plt.xlabel('Время')
plt.ylabel('Амплитуда (дБ)')
plt.title('Интерполированный с помощью таблицы сигнал сигнал sinc в дБ')
plt.legend()
plt.grid(True)
plt.ylim(-50, 5)
plt.show()

#Нахождение точек пересечения уровня половинной мощности и интерполированного сигнала
def razreshenie(sinc_interp, t_interp):
    half_P = -3
    crossings = []
    for i in range(1, len(sinc_interp)):
        if (sinc_interp[i - 1] - half_P) * (sinc_interp[i] - half_P) <= 0:
            t_cross = t_interp[i - 1] + (t_interp[i] - t_interp[i - 1]) * (half_P - sinc_interp[i - 1]) / (
                        sinc_interp[i] - sinc_interp[i - 1])
            crossings.append(t_cross)
    left_crossings = [x for x in crossings if x < 0]
    right_crossings = [x for x in crossings if x > 0]
    wl = (max(left_crossings)+min(left_crossings))/2
    wr = (min(right_crossings)+max(right_crossings))/2
    w = wr - wl
    return w

w1 = razreshenie(sinc_interp_1, t_interp_1)
w2 = razreshenie(sinc_interp_2, t_interp_2)
print(w1, w2)