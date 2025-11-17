der+1]
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
    print("5. Линейная интерполяция - самая быстрая, но наименее точная")import numpy as np
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

def calculate_mse(t_original, signal_original_db, t_interp, signal_interp_db):
    """Вычисление среднеквадратичной ошибки между интерполированным и теоретическим сигналом"""
    
    # Создаем теоретический sinc сигнал на интерполированной сетке
    t_theoretical = np.linspace(t_original[0], t_original[-1], len(t_interp))
    ideal_sinc = np.sinc(t_theoretical)
    ideal_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Вычисляем MSE в dB
    mse_db = np.mean((signal_interp_db - ideal_db) ** 2)
    
    return mse_db

def test_two_methods_comparison():
    """Сравнение двух методов (Wong и FFT) на чистом и зашумленном сигнале"""
    
    # Параметры сигнала
    t_original = np.linspace(-3, 3, 30)
    interp_factor = 4
    
    # Создаем чистый и зашумленный сигналы
    ideal_sinc = np.sinc(t_original)
    clean_db = 20 * np.log10(np.abs(ideal_sinc) + 1e-12)
    
    # Зашумленный сигнал
    np.random.seed(42)
    noise_level = 0.05
    noisy_sinc = ideal_sinc + np.random.normal(0, noise_level, len(ideal_sinc))
    noisy_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
    
    # Методы для сравнения
    methods = [
        ('Sinc (Wong)', sinc_interpolation_wong),
        ('FFT (без окна)', fourier_interpolation_basic)
    ]
    
    # Создаем графики
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Теоретический сигнал для сравнения
    t_dense = np.linspace(-3, 3, 1000)
    ideal_dense = np.sinc(t_dense)
    ideal_dense_db = 20 * np.log10(np.abs(ideal_dense) + 1e-12)
    
    # Цвета для методов
    colors = ['red', 'blue']
    
    # Тестируем на чистом сигнале
    print("ЧИСТЫЙ SINC-СИГНАЛ:")
    print("-" * 40)
    
    for idx, (method_name, method_func) in enumerate(methods):
        # Применяем интерполяцию к чистому сигналу
        t_interp, signal_interp = method_func(t_original, clean_db, interp_factor)
        
        # Вычисляем MSE
        mse_clean = calculate_mse(t_original, clean_db, t_interp, signal_interp)
        
        # Визуализация чистого сигнала
        axes[0, 0].plot(t_interp, signal_interp, 
                       color=colors[idx], 
                       linewidth=2, 
                       label=f'{method_name} (MSE: {mse_clean:.4f})',
                       alpha=0.8)
        
        print(f"{method_name:<15} | MSE: {mse_clean:.6f}")
    
    # Настройка графика для чистого сигнала
    axes[0, 0].plot(t_dense, ideal_dense_db, 'k-', linewidth=1, label='Теоретический sinc', alpha=0.5)
    axes[0, 0].plot(t_original, clean_db, 'ko', markersize=4, label='Исходные точки', alpha=0.8)
    axes[0, 0].set_title('Чистый sinc-сигнал')
    axes[0, 0].set_ylabel('Амплитуда (dB)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_ylim(-50, 5)
    
    # Тестируем на зашумленном сигнале
    print("\nЗАШУМЛЕННЫЙ SINC-СИГНАЛ:")
    print("-" * 40)
    
    for idx, (method_name, method_func) in enumerate(methods):
        # Применяем интерполяцию к зашумленному сигналу
        t_interp, signal_interp = method_func(t_original, noisy_db, interp_factor)
        
        # Вычисляем MSE относительно идеального sinc (без шума)
        mse_noisy = calculate_mse(t_original, clean_db, t_interp, signal_interp)
        
        # Визуализация зашумленного сигнала
        axes[0, 1].plot(t_interp, signal_interp, 
                       color=colors[idx], 
                       linewidth=2, 
                       label=f'{method_name} (MSE: {mse_noisy:.4f})',
                       alpha=0.8)
        
        print(f"{method_name:<15} | MSE: {mse_noisy:.6f}")
    
    # Настройка графика для зашумленного сигнала
    axes[0, 1].plot(t_dense, ideal_dense_db, 'k-', linewidth=1, label='Теоретический sinc', alpha=0.5)
    axes[0, 1].plot(t_original, noisy_db, 'ko', markersize=4, label='Исходные точки', alpha=0.8)
    axes[0, 1].set_title('Зашумленный sinc-сигнал')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_ylim(-50, 5)
    
    # Сравнение MSE на одном графике
    mse_clean_results = []
    mse_noisy_results = []
    method_names = []
    
    for method_name, method_func in methods:
        # Чистый сигнал
        t_interp_clean, signal_interp_clean = method_func(t_original, clean_db, interp_factor)
        mse_clean = calculate_mse(t_original, clean_db, t_interp_clean, signal_interp_clean)
        mse_clean_results.append(mse_clean)
        
        # Зашумленный сигнал  
        t_interp_noisy, signal_interp_noisy = method_func(t_original, noisy_db, interp_factor)
        mse_noisy = calculate_mse(t_original, clean_db, t_interp_noisy, signal_interp_noisy)
        mse_noisy_results.append(mse_noisy)
        
        method_names.append(method_name)
    
    # График сравнения MSE
    x = np.arange(len(methods))
    width = 0.35
    
    axes[1, 0].bar(x - width/2, mse_clean_results, width, label='Чистый сигнал', color='lightblue', alpha=0.8)
    axes[1, 0].bar(x + width/2, mse_noisy_results, width, label='Зашумленный сигнал', color='lightcoral', alpha=0.8)
    
    axes[1, 0].set_xlabel('Метод интерполяции')
    axes[1, 0].set_ylabel('MSE')
    axes[1, 0].set_title('Сравнение MSE для разных методов')
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(method_names, rotation=15)
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Разность MSE (шум - чистый)
    mse_difference = [noisy - clean for clean, noisy in zip(mse_clean_results, mse_noisy_results)]
    
    axes[1, 1].bar(method_names, mse_difference, color=['red', 'blue'], alpha=0.7)
    axes[1, 1].set_xlabel('Метод интерполяции')
    axes[1, 1].set_ylabel('ΔMSE (шум - чистый)')
    axes[1, 1].set_title('Влияние шума на точность интерполяции')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Добавляем значения на столбцы
    for i, v in enumerate(mse_difference):
        axes[1, 1].text(i, v + 0.001, f'{v:.4f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.show()
    
    # Вывод результатов
    print("\n" + "="*50)
    print("РЕЗЮМЕ:")
    print("="*50)
    
    for i, method_name in enumerate(method_names):
        print(f"{method_name}:")
        print(f"  MSE (чистый): {mse_clean_results[i]:.6f}")
        print(f"  MSE (шум):    {mse_noisy_results[i]:.6f}")
        print(f"  Ухудшение:    {mse_difference[i]:.6f}")
        print()
    
    # Определяем лучший метод для каждого случая
    best_clean_idx = np.argmin(mse_clean_results)
    best_noisy_idx = np.argmin(mse_noisy_results)
    
    print(f"Лучший метод для чистого сигнала: {method_names[best_clean_idx]}")
    print(f"Лучший метод для зашумленного сигнала: {method_names[best_noisy_idx]}")
    
    return {
        'methods': method_names,
        'mse_clean': mse_clean_results,
        'mse_noisy': mse_noisy_results,
        'mse_difference': mse_difference
    }

# Дополнительная функция для демонстрации артефактов
def show_artifacts_comparison():
    """Демонстрация артефактов интерполяции на зашумленном сигнале"""
    
    # Создаем зашумленный сигнал с более высоким уровнем шума
    t_original = np.linspace(-3, 3, 25)  # Меньше точек для наглядности
    ideal_sinc = np.sinc(t_original)
    
    np.random.seed(123)
    noise_level = 0.1  # Высокий уровень шума
    noisy_sinc = ideal_sinc + np.random.normal(0, noise_level, len(ideal_sinc))
    noisy_db = 20 * np.log10(np.abs(noisy_sinc) + 1e-12)
    
    # Интерполируем обоими методами
    t_wong, signal_wong = sinc_interpolation_wong(t_original, noisy_db, 6)  # Больший коэффициент
    t_fft, signal_fft = fourier_interpolation_basic(t_original, noisy_db, 6)
    
    # Теоретический сигнал
    t_dense = np.linspace(-3, 3, 1000)
    ideal_dense_db = 20 * np.log10(np.abs(np.sinc(t_dense)) + 1e-12)
    
    # Визуализация
    plt.figure(figsize=(12, 8))
    
    plt.plot(t_dense, ideal_dense_db, 'k-', linewidth=1, label='Теоретический sinc', alpha=0.5)
    plt.plot(t_original, noisy_db, 'ko', markersize=6, label='Зашумленные точки', alpha=0.8)
    
    plt.plot(t_wong, signal_wong, 'r-', linewidth=2, label='Sinc (Wong)', alpha=0.8)
    plt.plot(t_fft, signal_fft, 'b-', linewidth=2, label='FFT (без окна)', alpha=0.8)
    
    plt.xlabel('Время')
    plt.ylabel('Амплитуда (dB)')
    plt.title('Сравнение артефактов интерполяции на сильно зашумленном сигнале')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(-50, 5)
    plt.tight_layout()
    plt.show()
    
    # Вычисляем MSE для этого случая
    mse_wong = calculate_mse(t_original, noisy_db, t_wong, signal_wong)
    mse_fft = calculate_mse(t_original, noisy_db, t_fft, signal_fft)
    
    print("Сильно зашумленный сигнал (шум = 0.1):")
    print(f"Sinc (Wong) MSE: {mse_wong:.6f}")
    print(f"FFT (без окна) MSE: {mse_fft:.6f}")

# Запуск сравнения
if __name__ == "__main__":
    print("СРАВНЕНИЕ МЕТОДОВ ИНТЕРПОЛЯЦИИ: Sinc (Wong) vs FFT (без окна)")
    print("="*60)
    
    # Основное сравнение
    results = test_two_methods_comparison()
    
    print("\n" + "="*60)
    print("ДЕМОНСТРАЦИЯ АРТЕФАКТОВ ПРИ ВЫСОКОМ УРОВНЕ ШУМА")
    print("="*60)
    
    # Демонстрация артефактов
    show_artifacts_comparison()
    
    print("\n" + "="*60)
    print("ВЫВОДЫ:")
    print("="*60)
    print("1. Sinc (Wong) обычно дает лучшую точность на чистых сигналах")
    print("2. FFT метод более устойчив к шуму и дает меньше артефактов")
    print("3. При высоком уровне шума FFT метод часто предпочтительнее")
    print("4. Выбор метода зависит от уровня шума в ваших данных")
