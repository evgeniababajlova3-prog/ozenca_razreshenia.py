import matplotlib.pyplot as plt
import numpy as np


def plot_results(t_original, sinc_db, t_interp, sinc_interp, wl, wr, width, classical_pslr):
    """Упрощенная визуализация сигнала с основными параметрами."""
    plt.figure(figsize=(12, 6))

    # Исходный и интерполированный сигнал
    plt.plot(t_original, sinc_db, 'b-', linewidth=1, label='Исходный сигнал', alpha=0.7)
    plt.plot(t_interp, sinc_interp, 'g-', linewidth=1.5, label='Интерполированный сигнал', alpha=0.8)

    # Показываем ширину главного лепестка
    if wl is not None and wr is not None:
        plt.axvline(x=wl, color='red', linestyle='--', alpha=0.7, label=f'Ширина: {width:.3f}')
        plt.axvline(x=wr, color='red', linestyle='--', alpha=0.7)

        # Показываем максимальный УБЛ
        if classical_pslr > -80:
            pslr_level = -classical_pslr  # уровень бокового лепестка относительно главного
            plt.axhline(y=pslr_level, color='purple', linestyle=':',
                        alpha=0.7, label=f'УБЛ: {classical_pslr:.2f} дБ')

    plt.xlabel('Время/расстояние')
    plt.ylabel('Амплитуда (дБ)')
    plt.title('Анализ сигнала')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(-50, 5)
    plt.tight_layout()
    plt.show()