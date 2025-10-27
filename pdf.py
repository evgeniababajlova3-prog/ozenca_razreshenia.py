
import os
import subprocess
import datetime
import numpy as np
import matplotlib.pyplot as plt

# Сохраняем основные параметры РЛИ
with open(os.path.join(output_folder, "radar_params.txt"), 'w', encoding='utf-8') as f:
    f.write(f"Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей\n")
    f.write(f"Количество целей: {len(detected_peaks)}\n")
    f.write(f"Максимальная амплитуда РЛИ: {np.max(radar_image):.6f}\n")
    f.write(f"Минимальная амплитуда РЛИ: {np.min(radar_image):.6f}\n")
    f.write(f"Средняя амплитуда РЛИ: {np.mean(radar_image):.6f}\n")

# 1. Сохраняем основное РЛИ с целями
plt.figure(figsize=(12, 10))
plt.imshow(radar_image, cmap='hot', interpolation='nearest')
plt.colorbar(label='Амплитуда')
plt.title('Радиолокационное изображение с обнаруженными целями')
plt.xlabel('Пиксели по дальности')
plt.ylabel('Пиксели по азимуту')

# Отмечаем обнаруженные цели
for i, (y, x) in enumerate(detected_peaks):
    plt.plot(x, y, 'bx', markersize=10, markeredgewidth=2,
             label='Обнаруженные цели' if i == 0 else "")

plt.legend()
radar_image_path = os.path.join(output_folder, "radar_image.png")
plt.savefig(radar_image_path, dpi=150, bbox_inches='tight')
plt.close()

# 2. Сохраняем данные и графики для каждой цели
for i, target_data in enumerate(targets_data):
    target_id = i + 1
    target_coords = detected_peaks[i]

    # Создаем папку для цели
    target_folder = os.path.join(output_folder, f"target_{target_id}")
    os.makedirs(target_folder, exist_ok=True)

    # Сохраняем параметры цели
    with open(os.path.join(target_folder, f"target_{target_id}_params.txt"), 'w', encoding='utf-8') as f:
        f.write(f"Цель №{target_id}\n")
        f.write(f"Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}\n\n")

        f.write("Горизонтальное сечение:\n")
        f.write(f"Ширина главного лепестка: {target_data['h_width']:.4f}\n")
        f.write(f"Максимальный УБЛ: {target_data['h_pslr']:.2f} дБ\n")
        f.write(f"Интегральный УБЛ: {target_data['h_i_pslr']:.2f} дБ\n\n")

        f.write("Вертикальное сечение:\n")
        f.write(f"Ширина главного лепестка: {target_data['v_width']:.4f}\n")
        f.write(f"Максимальный УБЛ: {target_data['v_pslr']:.2f} дБ\n")
        f.write(f"Интегральный УБЛ: {target_data['v_i_pslr']:.2f} дБ\n")

    # Сохраняем массивы данных как в примере коллеги
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_window_linear.txt"),
               target_data['window_linear'], fmt='%.6f')
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_horizontal_section.txt"),
               np.column_stack([target_data['h_t'], target_data['h_signal_db']]), fmt='%.6f')
    np.savetxt(os.path.join(target_folder, f"target_{target_id}_vertical_section.txt"),
               np.column_stack([target_data['v_t'], target_data['v_signal_db']]), fmt='%.6f')

    # Окно в линейном масштабе
    plt.figure(figsize=(8, 6))
    plt.imshow(target_data['window_linear'], cmap='hot', interpolation='nearest')
    plt.colorbar(label='Амплитуда')
    plt.title(f'Цель {target_id} - Линейный масштаб')
    plt.xlabel('Отсчеты по дальности')
    plt.ylabel('Отсчеты по азимуту')
    linear_path = os.path.join(target_folder, f"target_{target_id}_linear.png")


plt.savefig(linear_path, dpi=150, bbox_inches='tight')
plt.close()

# Окно в логарифмическом масштабе
plt.figure(figsize=(8, 6))
window_db = 20 * np.log10(np.abs(target_data['window_linear']) + 1e-12)
plt.imshow(window_db, cmap='hot', interpolation='nearest')
plt.colorbar(label='Амплитуда (дБ)')
plt.title(f'Цель {target_id} - Логарифмический масштаб')
plt.xlabel('Отсчеты по дальности')
plt.ylabel('Отсчеты по азимуту')
db_path = os.path.join(target_folder, f"target_{target_id}_db.png")
plt.savefig(db_path, dpi=150, bbox_inches='tight')
plt.close()

# Горизонтальное сечение
plt.figure(figsize=(10, 5))
plt.plot(target_data['h_t'], target_data['h_signal_db'], 'b-', linewidth=2, alpha=0.8, label='Горизонтальное сечение')
if target_data['h_wl'] is not None and target_data['h_wr'] is not None:
    plt.axvline(x=target_data['h_wl'], color='red', linestyle='--', alpha=0.7, linewidth=2)
    plt.axvline(x=target_data['h_wr'], color='red', linestyle='--', alpha=0.7, linewidth=2)
    plt.axhline(y=-3, color='green', linestyle=':', alpha=0.7, linewidth=1.5)
    plt.text(0.05, 0.95, f'Ширина: {target_data["h_width"]:.3f}',
             transform=plt.gca().transAxes, fontsize=12,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
plt.xlabel('Время/расстояние')
plt.ylabel('Амплитуда (дБ)')
plt.title(f'Цель {target_id} - Горизонтальное сечение')
plt.grid(True, alpha=0.3)
plt.ylim(-50, 5)
plt.legend()
h_section_path = os.path.join(target_folder, f"target_{target_id}_horizontal.png")
plt.savefig(h_section_path, dpi=150, bbox_inches='tight')
plt.close()

# Вертикальное сечение
plt.figure(figsize=(10, 5))
plt.plot(target_data['v_t'], target_data['v_signal_db'], 'g-', linewidth=2, alpha=0.8, label='Вертикальное сечение')
if target_data['v_wl'] is not None and target_data['v_wr'] is not None:
    plt.axvline(x=target_data['v_wl'], color='red', linestyle='--', alpha=0.7, linewidth=2)
    plt.axvline(x=target_data['v_wr'], color='red', linestyle='--', alpha=0.7, linewidth=2)
    plt.axhline(y=-3, color='green', linestyle=':', alpha=0.7, linewidth=1.5)
    plt.text(0.05, 0.95, f'Ширина: {target_data["v_width"]:.3f}',
             transform=plt.gca().transAxes, fontsize=12,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
plt.xlabel('Время/расстояние')
plt.ylabel('Амплитуда (дБ)')
plt.title(f'Цель {target_id} - Вертикальное сечение')
plt.grid(True, alpha=0.3)
plt.ylim(-50, 5)
plt.legend()
v_section_path = os.path.join(target_folder, f"target_{target_id}_vertical.png")
plt.savefig(v_section_path, dpi=150, bbox_inches='tight')
plt.close()



def generate_typst_report(radar_image, detected_peaks, targets_data, output_folder):
    """
    Генерирует typst-код для отчета в точности как в примере коллеги.
    """
    current_time = datetime.datetime.now().strftime("%d.%m.%Y %H:%M")

    typ_content = f'''#set page(width: auto, height: auto, margin: 1.5cm)
#set text(font: "New Computer Modern", size: 12pt)
#show heading: set text(weight: "bold")

#align(center)[
  #text(size: 24pt, weight: "bold")[Анализ радиолокационного изображения]
  #text(size: 10pt)[  
  Отчет сгенерирован {current_time}
  ]
]

#block[
  Размер голограммы: {radar_image.shape[1]} × {radar_image.shape[0]} пикселей \
  Количество целей: {len(detected_peaks)}
]

// Визуализация РЛИ
#figure(
  image("{os.path.join(output_folder, "radar_image.png")}", width: 100%),
  caption: [Радиолокационное изображение с обнаруженными целями]
)

#pagebreak()
'''


# Добавляем данные по каждой цели как в примере
for i, target_data in enumerate(targets_data):
    target_id = i + 1
    target_coords = detected_peaks[i]
    target_folder = os.path.join(output_folder, f"target_{target_id}")

    typ_content += f'''
#align(center)[
  #text(size: 18pt, weight: "bold")[Цель №{target_id}]
]

Координаты: азимут {target_coords[0]}, дальность {target_coords[1]}

// Визуализация окна цели
#figure(
  grid(
    columns: 2,
    gutter: 1cm,
    [
      #figure(
        image("{os.path.join(target_folder, f"target_{target_id}_linear.png")}", width: 100%),
        caption: [Окно цели - линейный масштаб]
      )
    ],
    [
      #figure(
        image("{os.path.join(target_folder, f"target_{target_id}_db.png")}", width: 100%),
        caption: [Окно цели - логарифмический масштаб]
      )
    ]
  )
)

// Сечения цели
#figure(
  grid(
    columns: 2,
    gutter: 1cm,
    [
      #figure(
        image("{os.path.join(target_folder, f"target_{target_id}_horizontal.png")}", width: 100%),
        caption: [Горизонтальное сечение]
      )
    ],
    [
      #figure(
        image("{os.path.join(target_folder, f"target_{target_id}_vertical.png")}", width: 100%),
        caption: [Вертикальное сечение]
      )
    ]
  )
)

// Таблица параметров отклика от цели №{target_id}
#table(
  columns: 2,
  align: center,
  stroke: (x: 0.5pt, y: 0.5pt),
  inset: 5pt,

  [*Параметр*], [*Значение*],

  [Ширина главного лепестка (гориз.)], [{target_data['h_width']:.4f}],
  [Максимальный УБЛ (гориз.)], [{target_data['h_pslr']:.2f} дБ],
  [Интегральный УБЛ (гориз.)], [{target_data['h_i_pslr']:.2f} дБ],

  [Ширина главного лепестка (верт.)], [{target_data['v_width']:.4f}],
  [Максимальный УБЛ (верт.)], [{target_data['v_pslr']:.2f} дБ],
  [Интегральный УБЛ (верт.)], [{target_data['v_i_pslr']:.2f} дБ],
)

#pagebreak()
'''

# Сохраняем typst-файл
typ_filename = os.path.join(output_folder, "report.typ")
with open(typ_filename, 'w', encoding='utf-8') as f:
    f.write(typ_content)

return typ_filename


def compile_typst_to_pdf(typ_filename):
    """
    Компилирует typst-файл в PDF как в примере коллеги.
    """
    # Путь к компилятору typst
    typst_path = r'C:\Users\Gisich_AV\Desktop\typst\typst.exe'

    pdf_path = typ_filename.replace('.typ', '.pdf')

    try:
        result = subprocess.run([typst_path, 'compile', typ_filename, pdf_path],
                                capture_output=True, text=True, check=True)
        print(f"✓ PDF успешно сгенерирован: {pdf_path}")
        return pdf_path
    except subprocess.CalledProcessError as e:
        print(f"✗ Ошибка компиляции Typst: {e}")
        print(f"Stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print("✗ Компилятор Typst не найден. Убедитесь, что путь указан правильно.")
        return None


def clean_temp_files(output_folder):
    """
    Удаляет временные файлы как в примере коллеги.
    """
    # Сохраняем только PDF и основные изображения, остальное удаляем
    files_to_keep = ['report.pdf', 'radar_image.png']

    for root, dirs, files in os.walk(output_folder):
        for file in files:
            file_path = os.path.join(root, file)
            if file not in files_to_keep and not file.endswith('.pdf'):
                os.remove(file_path)

    print("✓ Временные файлы удалены")


def generate_analysis_report(radar_image, detected_peaks, targets_data, output_folder="analysis_report"):
    """
    Основная функция для генерации отчета в стиле примера коллеги.
    """
    print("Формирование отчета...")

    # Сохраняем все данные и графики
    plots_folder = save_report_data(radar_image, detected_peaks, targets_data, output_folder)

    # Генерируем typst-код



typ_filename = generate_typst_report(radar_image, detected_peaks, targets_data, plots_folder)

# Компилируем в PDF
pdf_path = compile_typst_to_pdf(typ_filename)

if pdf_path:
    print(f"✓ Отчет успешно создан: {pdf_path}")
    # Удаляем временные файлы
    clean_temp_files(output_folder)
else:
    print("✗ Не удалось создать отчет")

return pdf_path