import numpy as np

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