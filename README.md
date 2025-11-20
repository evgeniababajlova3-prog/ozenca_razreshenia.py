import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import firwin, lfilter
import warnings
warnings.filterwarnings('ignore')

def sinc_interpolation_wong(t_original, signal_db, interp_factor=4):
    """Sinc-интерполяция по алгоритму Wong 2005"""
    
    kernel_size = 8
    subsamples = 16
    
    def generate_coefficients(kernel_size, subsamples):
        coefficients = np.zeros((subsamples, kernel_size))
        window = np.kaiser(kernel_size, beta=2.5)

        for phase in range(subsamples):
            shift = phase / subsamples
            for k in range(kernel_size):
                pos = (k - (kernel_size-1)//2) - shift
                if abs(pos) < 1e-12:
                    sinc_val = 1.0
                else:
                    sinc_val = np.sin(np.pi * pos) / (np.pi * pos)
                coefficients[phase, k] = sinc_val * window[k]

            sum_coeff = np.sum(coefficients[phase])
            if abs(sum_coeff) > 1e-12:
                coefficients[phase] /= sum_coeff
        return coefficients

    # Основная функция интерполяции
    coefficients_table = generate_coefficients(kernel_size, subsamples)
    signal_linear = 10 ** (signal_db / 20)
    dt = t_original[1] - t_original[0]
    
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    signal_interp = np.zeros_like(t_interp)
    half_kernel = kernel_size // 2

    for i, t_point in enumerate(t_interp):
        idx_center = int(np.floor((t_point - t_original[0]) / dt))
        fractional = (t_point - t_original[idx_center]) / dt
        
        phase_idx = int(fractional * subsamples)
        phase_idx = max(0, min(phase_idx, subsamples - 1))
        coeffs = coefficients_table[phase_idx]

        start_idx = idx_center - half_kernel
        end_idx = start_idx + kernel_size

        if start_idx < 0 or end_idx > len(signal_linear):
            signal_interp[i] = np.interp(t_point, t_original, signal_linear)
        else:
            points = signal_linear[start_idx:end_idx]
            signal_interp[i] = np.sum(points * coeffs)

    signal_interp_db = 20 * np.log10(np.abs(signal_interp) + 1e-12)
    return t_interp, signal_interp_db

def fourier_interpolation_basic(t_original, signal_db, interp_factor=4):
    """Базовая Фурье-интерполяция с дополнением нулями"""
    
    signal_linear = 10 ** (signal_db / 20)
    N_original = len(signal_linear)
    
    spectrum = np.fft.fft(signal_linear)
    N_new = N_original * interp_factor
    new_spectrum = np.zeros(N_new, dtype=complex)
    
    half_original = N_original // 2
    new_spectrum[:half_original] = spectrum[:half_original]
    new_spectrum[-half_original:] = spectrum[-half_original:]
    new_spectrum *= interp_factor
    
    signal_interp_linear = np.fft.ifft(new_spectrum).real
    t_interp = np.linspace(t_original[0], t_original[-1], N_new)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def fourier_interpolation_hann(t_original, signal_db, interp_factor=4):
    """Фурье-интерполяция с оконной функцией Hann"""
    
    signal_linear = 10 ** (signal_db / 20)
    N_original = len(signal_linear)
    
    # Оконная функция
    window = np.hanning(N_original)
    signal_windowed = signal_linear * window
    
    spectrum = np.fft.fft(signal_windowed)
    N_new = N_original * interp_factor
    new_spectrum = np.zeros(N_new, dtype=complex)
    
    half_original = N_original // 2
    new_spectrum[:half_original] = spectrum[:half_original]
    new_spectrum[-half_original:] = spectrum[-half_original:]
    new_spectrum *= interp_factor
    
    signal_interp_linear = np.fft.ifft(new_spectrum).real
    
    # Компенсация оконной функции
    window_interp = np.interp(
        np.linspace(0, N_original-1, N_new),
        np.arange(N_original),
        window
    )
    
    signal_interp_linear = np.divide(
        signal_interp_linear, 
        window_interp + 1e-12,
        out=signal_interp_linear,
        where=window_interp > 1e-6
    )
    
    t_interp = np.linspace(t_original[0], t_original[-1], N_new)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def cubic_spline_interpolation(t_original, signal_db, interp_factor=4):
    """Кубическая сплайн-интерполяция"""
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    
    cs = CubicSpline(t_original, signal_linear)
    signal_interp_linear = cs(t_interp)
    
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    return t_interp, signal_interp_db

def linear_interpolation(t_original, signal_db, interp_factor=4):
    """Линейная интерполяция"""
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    
    signal_interp_linear = np.interp(t_interp, t_original, signal_linear)
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def lagrange_interpolation(t_original, signal_db, interp_factor=4, order=3):
    """Полиномиальная интерполяция методом Лагранжа"""
    
    def lagrange_polynomial(x, x_points, y_points):
        total = 0
        n = len(x_points)
        for i in range(n):
            xi, yi = x_points[i], y_points[i]
            product = 1
            for j in range(n):
                if i != j:
                    xj = x_points[j]
                    product *= (x - xj) / (xi - xj)
            total += yi * product
        return total
    
    signal_linear = 10 ** (signal_db / 20)
    t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
    signal_interp_linear = np.zeros_like(t_interp)
    
    # Применяем скользящее окно Лагранжа
    for i, t_point in enumerate(t_interp):
        # Находим ближайшие точки
        distances = np.abs(t_original - t_point)
        nearest_indices = np.argsort(distances)[:order+1]
        nearest_t = t_original[nearest_indices]
        nearest_y = signal_linear[nearest_indices]
        
        signal_interp_linear[i] = lagrange_polynomial(t_point, nearest_t, nearest_y)
    
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    return t_interp, signal_interp_db

def fir_filter_interpolation(t_original, signal_db, interp_factor=4):
    """Интерполяция через FIR фильтр"""
    
    signal_linear = 10 ** (signal_db / 20)
    
    # Создаем FIR фильтр для интерполяции
    numtaps = 31
    cutoff = 0.8 / interp_factor
    fir_coeff = firwin(numtaps, cutoff)
    
    # Увеличиваем частоту дискретизации
    upsampled = np.zeros(len(signal_linear) * interp_factor)
    upsampled[::interp_factor] = signal_linear
    
    # Применяем FIR фильтр
    signal_interp_linear = lfilter(fir_coeff, 1.0, upsampled)
    
    # Обрезаем переходный процесс
    signal_interp_linear = signal_interp_linear[numtaps//2:]
    
    t_interp = np.linspace(t_original[0], t_original[-1], len(signal_interp_linear))
    signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
    
    return t_interp, signal_interp_db

def calculate_mse(t_original, signal_original_db, t_interp, signal_interp_db):
    """Вычисление среднеквадратичной ошибки между интерполированным и теоретическим сигналом"""
    
    # Создаем теоретический sinc сигнал на интерполированной сетке
    t_theoretical = np.linspace(t_original[0], t_original[-1], len(t_interp))
    ideal_sinc = np.sinc(t_theoretical)
    ideal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Вычисляем MSE в dB
    mse_db = np.mean((signal_interp_db - ideal_db) ** 2)
    
    # Также вычисляем MSE в линейной области для сравнения
    signal_interp_linear = 10 ** (signal_interp_db / 20)
    ideal_linear = np.abs(np.sinc(t_theoretical))
    mse_linear = np.mean((signal_interp_linear - ideal_linear) ** 2)
    
    return mse_db, mse_linear

def test_interpolation_methods():
    """Тестирование всех методов интерполяции на идеальном sinc-сигнале"""
    
    # Создаем идеальный sinc-сигнал для тестирования
    t_original = np.linspace(-3, 3, 30)  # 30 точек
    ideal_sinc = np.sinc(t_original)
    signal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Список методов для тестирования
    methods = [
        ('Sinc (Wong)', sinc_interpolation_wong),
        ('FFT (базовый)', fourier_interpolation_basic),
        ('FFT + Hann', fourier_interpolation_hann),
        ('Кубический сплайн', cubic_spline_interpolation),
        ('Линейная', linear_interpolation),
        ('Лагранж (порядок 3)', lambda t, s: lagrange_interpolation(t, s, 4, 3)),
        ('FIR фильтр', fir_filter_interpolation)
    ]
    
    # Результаты
    results = []
    
    # Создаем график для сравнения
    plt.figure(figsize=(15, 10))
    
    # Рисуем теоретический sinc
    t_dense = np.linspace(-3, 3, 1000)
    ideal_dense = np.sinc(t_dense)
    ideal_dense_db = 20 * np.log10(np.abs(ideal_dense) + 1e-12)
    plt.plot(t_dense, ideal_dense_db, 'k-', linewidth=3, label='Теоретический sinc', alpha=0.7)
    
    # Рисуем исходные точки
    plt.plot(t_original, signal_db, 'ko', markersize=8, label='Исходные точки', alpha=0.8)
    
    # Тестируем каждый метод
    colors = plt.cm.tab10(np.linspace(0, 1, len(methods)))
    
    for idx, (method_name, method_func) in enumerate(methods):
        try:
            # Применяем интерполяцию
            t_interp, signal_interp = method_func(t_original, signal_db, 4)
            
            # Вычисляем MSE
            mse_db, mse_linear = calculate_mse(t_original, signal_db, t_interp, signal_interp)
            
            # Сохраняем результаты
            results.append({
                'method': method_name,
                'mse_db': mse_db,
                'mse_linear': mse_linear,
                't_interp': t_interp,
                'signal_interp': signal_interp
            })
            
            # Рисуем интерполированный сигнал
            plt.plot(t_interp, signal_interp, 
                    color=colors[idx], 
                    linewidth=2, 
                    label=f'{method_name} (MSE: {mse_db:.4f})',
                    alpha=0.8)
            
            print(f"{method_name:<20} | MSE (dB): {mse_db:.6f} | MSE (linear): {mse_linear:.6f}")
            
        except Exception as e:
            print(f"{method_name:<20} | ОШИБКА: {str(e)}")
            results.append({
                'method': method_name,
                'mse_db': float('inf'),
                'mse_linear': float('inf'),
                'error': str(e)
            })
    
    # Настройка графика
    plt.xlabel('Время')
    plt.ylabel('Амплитуда (dB)')
    plt.title('Сравнение методов интерполяции на sinc-сигнале')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.ylim(-50, 5)
    plt.tight_layout()
    plt.show()
    
    # Сортируем методы по качеству (по MSE в dB)
    results_sorted = sorted([r for r in results if 'mse_db' in r], key=lambda x: x['mse_db'])
    
    # Выводим таблицу результатов
    print("\n" + "="*70)
    print("РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ МЕТОДОВ ИНТЕРПОЛЯЦИИ")
    print("="*70)
    print(f"{'Метод':<25} {'MSE (dB)':<15} {'MSE (linear)':<15}")
    print("-"*70)
    
    for result in results_sorted:
        print(f"{result['method']:<25} {result['mse_db']:<15.6f} {result['mse_linear']:<15.6f}")
    
    # Выводим лучший метод
    if results_sorted:
        best_method = results_sorted[0]
        print("-"*70)
        print(f"ЛУЧШИЙ МЕТОД: {best_method['method']}")
        print(f"MSE (dB): {best_method['mse_db']:.6f}")
        print(f"MSE (linear): {best_method['mse_linear']:.6f}")
    
    return results_sorted

def test_with_noisy_signal():
    """Дополнительный тест с зашумленным сигналом"""
    
    print("\n" + "="*70)
    print("ТЕСТ С ЗАШУМЛЕННЫМ SINC-СИГНАЛОМ")
    print("="*70)
    
    # Создаем зашумленный sinc-сигнал
    t_original = np.linspace(-3, 3, 30)
    ideal_sinc = np.sinc(t_original)
    
    # Добавляем шум
    np.random.seed(42)
    noise = np.random.normal(0, 0.05, len(ideal_sinc))
    noisy_sinc = ideal_sinc + noise
    signal_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
    
    # Тестируем только лучшие методы
    methods = [
        ('Sinc (Wong)', sinc_interpolation_wong),
        ('FFT + Hann', fourier_interpolation_hann),
        ('Кубический сплайн', cubic_spline_interpolation),
        ('FIR фильтр', fir_filter_interpolation)
    ]
    
    results = []
    
    for method_name, method_func in methods:
        try:
            t_interp, signal_interp = method_func(t_original, signal_db, 4)
            
            # Вычисляем MSE относительно идеального sinc (без шума)
            t_theoretical = np.linspace(t_original[0], t_original[-1], len(t_interp))
            ideal_sinc_dense = np.sinc(t_theoretical)
            ideal_db = 20 * np.log10(np.abs(ideal_sinc_dense) + 1e-12)
            
            mse_db = np.mean((signal_interp - ideal_db) ** 2)
            
            signal_interp_linear = 10 ** (signal_interp / 20)
            ideal_linear = np.abs(np.sinc(t_theoretical))
            mse_linear = np.mean((signal_interp_linear - ideal_linear) ** 2)
            
            results.append({
                'method': method_name,
                'mse_db': mse_db,
                'mse_linear': mse_linear
            })
            
            print(f"{method_name:<20} | MSE (dB): {mse_db:.6f} | MSE (linear): {mse_linear:.6f}")
            
        except Exception as e:
            print(f"{method_name:<20} | ОШИБКА: {str(e)}")
    
    # Сортируем по качеству
    results_sorted = sorted(results, key=lambda x: x['mse_db'])
    
    print("\nЛучший метод для зашумленного сигнала:")
    print(f"{results_sorted[0]['method']} с MSE (dB) = {results_sorted[0]['mse_db']:.6f}")
    
    return results_sorted

def compare_interpolation_methods(t_original, signal_db, method='auto'):
    """
    Универсальная функция для выбора метода интерполяции
    
    Parameters:
    -----------
    t_original : array
        Исходная временная ось
    signal_db : array
        Сигнал в dB
    method : str
        Метод интерполяции или 'auto' для автоматического выбора
    
    Returns:
    --------
    tuple : (t_interp, signal_interp_db, method_used)
    """
    
    if method == 'auto':
        # Простой автоматический выбор на основе тестов
        # В реальном коде можно добавить более сложную логику
        method = 'fft_hann'
    
    method_functions = {
        'sinc_wong': sinc_interpolation_wong,
        'fft_basic': fourier_interpolation_basic,
        'fft_hann': fourier_interpolation_hann,
        'cubic_spline': cubic_spline_interpolation,
        'linear': linear_interpolation,
        'lagrange': lambda t, s: lagrange_interpolation(t, s, 4, 3),
        'fir_filter': fir_filter_interpolation
    }
    
    if method not in method_functions:
        raise ValueError(f"Неизвестный метод: {method}")
    
    t_interp, signal_interp_db = method_functions[method](t_original, signal_db, 4)
    
    return t_interp, signal_interp_db, method

# Запуск тестирования
if __name__ == "__main__":
    print("ТЕСТИРОВАНИЕ МЕТОДОВ ИНТЕРПОЛЯЦИИ")
    print("Сравнение среднеквадратичного отклонения от теоретического sinc-сигнала")
    print("="*70)
    
    # Основной тест на идеальном sinc
    results = test_interpolation_methods()
    
    # Дополнительный тест с шумом
    noisy_results = test_with_noisy_signal()
    
    print("\n" + "="*70)
    print("ВЫВОДЫ:")
    print("="*70)
    print("1. FFT методы обычно показывают лучшую устойчивость к артефактам")
    print("2. Sinc (Wong) дает максимальную точность на идеальных сигналах") 
    print("3. FIR фильтр хорошо работает с зашумленными сигналами")
    print("4. Для большинства практических задач рекомендуется FFT + Hann")
    print("5. Линейная интерполяция - самая быстрая, но наименее точная")        
        # Идеальный sinc
        t_ideal = np.linspace(-3, 3, 30)
        ideal_sinc = np.sinc(t_ideal)
        ideal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
        
        # Зашумленный sinc
        np.random.seed(42)
        noise = np.random.normal(0, 0.05, len(ideal_sinc))
        noisy_sinc = ideal_sinc + noise
        noisy_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
        
        # Два близких пика (для проверки разрешения)
        t_peaks = np.linspace(-2, 2, 40)
        two_peaks = np.sinc(t_peaks * 2) + 0.7 * np.sinc((t_peaks - 0.5) * 2)
        peaks_db = 20 * np.log10(np.abs(two_peaks) + 1e-12)
        
        return {
            'ideal_sinc': (t_ideal, ideal_db),
            'noisy_sinc': (t_ideal, noisy_db),
            'two_peaks': (t_peaks, peaks_db)
        }
    
    def sinc_interpolation_wong(self, t_original, signal_db, interp_factor=4):
        """Sinc-интерполяция по алгоритму Wong 2005"""
        
        kernel_size = 8
        subsamples = 16
        
        def generate_coefficients(kernel_size, subsamples):
            coefficients = np.zeros((subsamples, kernel_size))
            window = np.kaiser(kernel_size, beta=2.5)

            for phase in range(subsamples):
                shift = phase / subsamples
                for k in range(kernel_size):
                    pos = (k - (kernel_size-1)//2) - shift
                    if abs(pos) < 1e-12:
                        sinc_val = 1.0
                    else:
                        sinc_val = np.sin(np.pi * pos) / (np.pi * pos)
                    coefficients[phase, k] = sinc_val * window[k]

                sum_coeff = np.sum(coefficients[phase])
                if abs(sum_coeff) > 1e-12:
                    coefficients[phase] /= sum_coeff
            return coefficients

        # Основная функция интерполяции
        coefficients_table = generate_coefficients(kernel_size, subsamples)
        signal_linear = 10 ** (signal_db / 20)
        dt = t_original[1] - t_original[0]
        
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        signal_interp = np.zeros_like(t_interp)
        half_kernel = kernel_size // 2

        for i, t_point in enumerate(t_interp):
            idx_center = int(np.floor((t_point - t_original[0]) / dt))
            fractional = (t_point - t_original[idx_center]) / dt
            
            phase_idx = int(fractional * subsamples)
            phase_idx = max(0, min(phase_idx, subsamples - 1))
            coeffs = coefficients_table[phase_idx]

            start_idx = idx_center - half_kernel
            end_idx = start_idx + kernel_size

            if start_idx < 0 or end_idx > len(signal_linear):
                signal_interp[i] = np.interp(t_point, t_original, signal_linear)
            else:
                points = signal_linear[start_idx:end_idx]
                signal_interp[i] = np.sum(points * coeffs)

        signal_interp_db = 20 * np.log10(np.abs(signal_interp) + 1e-12)
        return t_interp, signal_interp_db
    
    def fourier_interpolation_basic(self, t_original, signal_db, interp_factor=4):
        """Базовая Фурье-интерполяция с дополнением нулями"""
        
        signal_linear = 10 ** (signal_db / 20)
        N_original = len(signal_linear)
        
        spectrum = np.fft.fft(signal_linear)
        N_new = N_original * interp_factor
        new_spectrum = np.zeros(N_new, dtype=complex)
        
        half_original = N_original // 2
        new_spectrum[:half_original] = spectrum[:half_original]
        new_spectrum[-half_original:] = spectrum[-half_original:]
        new_spectrum *= interp_factor
        
        signal_interp_linear = np.fft.ifft(new_spectrum).real
        t_interp = np.linspace(t_original[0], t_original[-1], N_new)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def fourier_interpolation_hann(self, t_original, signal_db, interp_factor=4):
        """Фурье-интерполяция с оконной функцией Hann"""
        
        signal_linear = 10 ** (signal_db / 20)
        N_original = len(signal_linear)
        
        # Оконная функция
        window = np.hanning(N_original)
        signal_windowed = signal_linear * window
        
        spectrum = np.fft.fft(signal_windowed)
        N_new = N_original * interp_factor
        new_spectrum = np.zeros(N_new, dtype=complex)
        
        half_original = N_original // 2
        new_spectrum[:half_original] = spectrum[:half_original]
        new_spectrum[-half_original:] = spectrum[-half_original:]
        new_spectrum *= interp_factor
        
        signal_interp_linear = np.fft.ifft(new_spectrum).real
        
        # Компенсация оконной функции
        window_interp = np.interp(
            np.linspace(0, N_original-1, N_new),
            np.arange(N_original),
            window
        )
        
        signal_interp_linear = np.divide(
            signal_interp_linear, 
            window_interp + 1e-12,
            out=signal_interp_linear,
            where=window_interp > 1e-6
        )
        
        t_interp = np.linspace(t_original[0], t_original[-1], N_new)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def cubic_spline_interpolation(self, t_original, signal_db, interp_factor=4):
        """Кубическая сплайн-интерполяция"""
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        
        cs = CubicSpline(t_original, signal_linear)
        signal_interp_linear = cs(t_interp)
        
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        return t_interp, signal_interp_db
    
    def linear_interpolation(self, t_original, signal_db, interp_factor=4):
        """Линейная интерполяция"""
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        
        signal_interp_linear = np.interp(t_interp, t_original, signal_linear)
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def lagrange_interpolation(self, t_original, signal_db, interp_factor=4, order=3):
        """Полиномиальная интерполяция методом Лагранжа"""
        
        def lagrange_polynomial(x, x_points, y_points):
            total = 0
            n = len(x_points)
            for i in range(n):
                xi, yi = x_points[i], y_points[i]
                product = 1
                for j in range(n):
                    if i != j:
                        xj = x_points[j]
                        product *= (x - xj) / (xi - xj)
                total += yi * product
            return total
        
        signal_linear = 10 ** (signal_db / 20)
        t_interp = np.linspace(t_original[0], t_original[-1], len(t_original) * interp_factor)
        signal_interp_linear = np.zeros_like(t_interp)
        
        # Применяем скользящее окно Лагранжа
        for i, t_point in enumerate(t_interp):
            # Находим ближайшие точки
            distances = np.abs(t_original - t_point)
            nearest_indices = np.argsort(distances)[:order+1]
            nearest_t = t_original[nearest_indices]
            nearest_y = signal_linear[nearest_indices]
            
            signal_interp_linear[i] = lagrange_polynomial(t_point, nearest_t, nearest_y)
        
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        return t_interp, signal_interp_db
    
    def fir_filter_interpolation(self, t_original, signal_db, interp_factor=4):
        """Интерполяция через FIR фильтр (метод увеличения частоты дискретизации)"""
        
        signal_linear = 10 ** (signal_db / 20)
        
        # Создаем FIR фильтр для интерполяции
        numtaps = 31
        cutoff = 0.8 / interp_factor
        fir_coeff = firwin(numtaps, cutoff)
        
        # Увеличиваем частоту дискретизации
        upsampled = np.zeros(len(signal_linear) * interp_factor)
        upsampled[::interp_factor] = signal_linear
        
        # Применяем FIR фильтр
        signal_interp_linear = lfilter(fir_coeff, 1.0, upsampled)
        
        # Обрезаем переходный процесс
        signal_interp_linear = signal_interp_linear[numtaps//2:]
        
        t_interp = np.linspace(t_original[0], t_original[-1], len(signal_interp_linear))
        signal_interp_db = 20 * np.log10(np.abs(signal_interp_linear) + 1e-12)
        
        return t_interp, signal_interp_db
    
    def calculate_metrics(self, t_original, signal_original_db, t_interp, signal_interp_db, method_name):
        """Вычисление метрик качества интерполяции"""
        
        # Интерполируем исходные точки для сравнения
        signal_linear_original = 10 ** (signal_original_db / 20)
        signal_linear_interp = 10 ** (signal_interp_db / 20)
        
        # Находим соответствующие точки в интерполированном сигнале
        original_points_interp = np.interp(t_original, t_interp, signal_linear_interp)
        
        # Метрики
        mse = np.mean((signal_linear_original - original_points_interp) ** 2)
        
        # Плавность (средняя производная)
        derivative = np.diff(signal_interp_db)
        smoothness = np.mean(np.abs(derivative))
        
        # Сохранение энергии
        energy_original = np.sum(signal_linear_original ** 2)
        energy_interp = np.sum(signal_linear_interp ** 2) / len(signal_linear_interp) * len(signal_linear_original)
        energy_preservation = abs(energy_original - energy_interp) / energy_original
        
        # Максимальная ошибка
        max_error = np.max(np.abs(signal_linear_original - original_points_interp))
        
        return {
            'mse': mse,
            'smoothness': smoothness,
            'energy_preservation': energy_preservation,
            'max_error': max_error
        }
    
    def test_all_methods(self, t_original, signal_db, signal_name="Сигнал"):
        """Тестирование всех методов на заданном сигнале"""
        
        interp_factor = 4
        metrics = {}
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        axes = axes.flatten()
        
        methods = [
            ('sinc_wong', self.sinc_interpolation_wong),
            ('fft_basic', self.fourier_interpolation_basic),
            ('fft_hann', self.fourier_interpolation_hann),
            ('cubic_spline', self.cubic_spline_interpolation),
            ('linear', self.linear_interpolation),
            ('lagrange', lambda t, s: self.lagrange_interpolation(t, s, interp_factor, 3)),
            ('fir_filter', self.fir_filter_interpolation)
        ]
        
        for idx, (method_key, method_func) in enumerate(methods):
            if idx >= len(axes):
                break
                
            try:
                t_interp, signal_interp = method_func(t_original, signal_db, interp_factor)
                
                # Вычисляем метрики
                metrics[method_key] = self.calculate_metrics(
                    t_original, signal_db, t_interp, signal_interp, method_key
                )
                
                # Визуализация
                ax = axes[idx]
                ax.plot(t_original, signal_db, 'bo-', label='Исходный', markersize=4, linewidth=1, alpha=0.7)
                ax.plot(t_interp, signal_interp, 'r-', label='Интерполированный', linewidth=1.5, alpha=0.8)
                ax.set_title(f'{self.methods[method_key]}', fontsize=12)
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                
                # Добавляем метрики на график
                metric_text = (f'MSE: {metrics[method_key]["mse"]:.2e}\n'
                             f'Плавность: {metrics[method_key]["smoothness"]:.3f}\n'
                             f'Энергия: {metrics[method_key]["energy_preservation"]:.2e}')
                ax.text(0.02, 0.98, metric_text, transform=ax.transAxes, fontsize=8,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
                
            except Exception as e:
                print(f"Ошибка в методе {method_key}: {e}")
                axes[idx].text(0.5, 0.5, f'Ошибка:\n{str(e)}', 
                              transform=axes[idx].transAxes, ha='center', va='center')
                axes[idx].set_title(f'{self.methods[method_key]} - ОШИБКА')
        
        # Скрываем пустые subplots
        for idx in range(len(methods), len(axes)):
            axes[idx].set_visible(False)
        
        plt.suptitle(f'Сравнение методов интерполяции: {signal_name}', fontsize=16)
        plt.tight_layout()
        plt.show()
        
        return metrics
    
    def comprehensive_test(self):
        """Комплексное тестирование на всех типах сигналов"""
        
        test_signals = self.generate_test_signals()
        all_metrics = {}
        
        for signal_name, (t, signal_db) in test_signals.items():
            print(f"\n{'='*50}")
            print(f"Тестирование: {signal_name}")
            print(f"{'='*50}")
            
            metrics = self.test_all_methods(t, signal_db, signal_name)
            all_metrics[signal_name] = metrics
            
            # Сводная таблица метрик
            self.print_metrics_table(metrics, signal_name)
        
        return all_metrics
    
    def print_metrics_table(self, metrics, signal_name):
        """Печать сводной таблицы метрик"""
        
        print(f"\nСводные метрики для {signal_name}:")
        print("-" * 80)
        print(f"{'Метод':<20} {'MSE':<12} {'Плавность':<12} {'Энергия':<12} {'Макс. ошибка':<12}")
        print("-" * 80)
        
        for method_key, metric in metrics.items():
            print(f"{self.methods[method_key]:<20} {metric['mse']:<12.2e} {metric['smoothness']:<12.3f} "
                  f"{metric['energy_preservation']:<12.2e} {metric['max_error']:<12.2e}")
    
    def find_best_method(self, metrics_dict, criterion='mse'):
        """Нахождение лучшего метода по заданному критерию"""
        
        best_methods = {}
        
        for signal_name, metrics in metrics_dict.items():
            best_score = float('inf')
            best_method = None
            
            for method_key, metric in metrics.items():
                score = metric[criterion]
                if score < best_score:
                    best_score = score
                    best_method = method_key
            
            best_methods[signal_name] = (best_method, best_score)
        
        print(f"\nЛучшие методы по критерию '{criterion}':")
        print("-" * 50)
        for signal_name, (method, score) in best_methods.items():
            print(f"{signal_name:<15} -> {self.methods[method]:<20} (значение: {score:.2e})")
        
        return best_methods

def universal_interpolation_analysis(t_original, signal_db, time_range, method='auto', interp_factor=4):
    """
    Универсальная функция анализа с автоматическим выбором или указанием метода интерполяции
    
    Parameters:
    -----------
    t_original : array
        Исходная временная ось
    signal_db : array
        Сигнал в dB
    time_range : float
        Диапазон времени (для совместимости)
    method : str
        Метод интерполяции: 'auto', 'sinc_wong', 'fft_basic', 'fft_hann', 
                           'cubic_spline', 'linear', 'lagrange', 'fir_filter'
    interp_factor : int
        Коэффициент интерполяции
    
    Returns:
    --------
    dict : Результаты анализа
    """
    
    comparator = InterpolationComparator()
    
    # Автоматический выбор метода
    if method == 'auto':
        # Простой тест для автоматического выбора
        signal_linear = 10 ** (signal_db / 20)
        noise_level = np.std(signal_linear) / np.mean(signal_linear)
        
        if noise_level > 0.1:  # Зашумленный сигнал
            method = 'fft_hann'
        elif len(t_original) < 20:  # Мало точек
            method = 'cubic_spline'
        else:  # Стандартный случай
            method = 'sinc_wong'
    
    # Применение выбранного метода
    method_functions = {
        'sinc_wong': comparator.sinc_interpolation_wong,
        'fft_basic': comparator.fourier_interpolation_basic,
        'fft_hann': comparator.fourier_interpolation_hann,
        'cubic_spline': comparator.cubic_spline_interpolation,
        'linear': comparator.linear_interpolation,
        'lagrange': lambda t, s: comparator.lagrange_interpolation(t, s, interp_factor, 3),
        'fir_filter': comparator.fir_filter_interpolation
    }
    
    if method not in method_functions:
        raise ValueError(f"Неизвестный метод: {method}. Доступные: {list(method_functions.keys())}")
    
    t_interp, sinc_interp = method_functions[method](t_original, signal_db, interp_factor)
    
    # Анализ характеристик (упрощенный)
    def find_main_lobe_width(t, signal_db):
        center_idx = np.argmax(signal_db)
        threshold = signal_db[center_idx] - 3
        
        # Поиск точек пересечения -3 dB
        left_idx = center_idx
        while left_idx > 0 and signal_db[left_idx] > threshold:
            left_idx -= 1
        
        right_idx = center_idx
        while right_idx < len(signal_db) - 1 and signal_db[right_idx] > threshold:
            right_idx += 1
        
        if left_idx > 0 and right_idx < len(signal_db) - 1:
            # Линейная интерполяция
            wl = t[left_idx] + (t[left_idx+1] - t[left_idx]) * (threshold - signal_db[left_idx]) / (signal_db[left_idx+1] - signal_db[left_idx])
            wr = t[right_idx-1] + (t[right_idx] - t[right_idx-1]) * (threshold - signal_db[right_idx-1]) / (signal_db[right_idx] - signal_db[right_idx-1])
            width = wr - wl
            return wl, wr, width
        return None, None, 0
    
    wl, wr, width = find_main_lobe_width(t_interp, sinc_interp)
    
    # Упрощенный расчет УБЛ
    def calculate_sidelobe_levels(t, signal_db):
        signal_linear = 10 ** (signal_db / 20)
        main_lobe_max = np.max(signal_linear)
        main_lobe_idx = np.argmax(signal_linear)
        
        # Боковые лепестки - все кроме области вокруг главного лепестка
        sidelobe_mask = np.ones(len(signal_linear), dtype=bool)
        sidelobe_mask[max(0, main_lobe_idx-5):min(len(signal_linear), main_lobe_idx+6)] = False
        
        if np.any(sidelobe_mask):
            max_sidelobe = np.max(signal_linear[sidelobe_mask])
            classical_pslr = 20 * np.log10(max_sidelobe / main_lobe_max)
            
            # Упрощенный интегральный УБЛ
            power_main = np.sum(signal_linear[~sidelobe_mask] ** 2)
            power_sidelobes = np.sum(signal_linear[sidelobe_mask] ** 2)
            integral_pslr = 10 * np.log10(power_sidelobes / power_main)
        else:
            classical_pslr = integral_pslr = -80
        
        return classical_pslr, integral_pslr
    
    classica
