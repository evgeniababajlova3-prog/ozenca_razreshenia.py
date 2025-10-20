ozenca_razreshenia.py
import os
import os.path
import subprocess

column =["1","2","3"]
headers =["один","два","три"]

result_folder = "data_20_10_2025"
os.makedirs(result_folder, exist_ok=True)

typ_filename = os.path.join(result_folder, "table.typ")

table_data = f"""
            #set page(width: auto, height: auto, margin: 1cm)
            #align(center)[#text(size: 28pt, weight: "bold")[Таблица]]
    
            #figure(
        table(
            columns: 3,
            align: center,
            stroke: (x: 0.5pt, y: 0.5pt),
        
            // Заголовки столбцов
            [#text(size: 22pt, weight: "bold")[{headers[0]}]], 
            [#text(size: 22pt, weight: "bold")[{headers[1]}]],
            [#text(size: 22pt, weight: "bold")[{headers[2]}]],
            
            // Строки данных
            [#text(size: 22pt, weight: "bold")[{column[0]}]],
            [#text(size: 22pt, weight: "bold")[{column[1]}]],
            [#text(size: 22pt, weight: "bold")[{column[2]}]],
                ),
            )
            """
with open(typ_filename, 'w', encoding='utf-8') as typ_file:
    typ_file.write(table_data)
typ_path = os.path.join(result_folder, "table.typ")
pdf_path = typ_path.replace('.typ', '.pdf')

            # Вызываем компилятор Typst
typst_path = r'C:\Users\Gisich_AV\Desktop\typst\typst\typst.exe'
result = subprocess.run([typst_path, 'compile', typ_path, pdf_path])

import os
import shutil
import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from nir_rlx_reader.reader import NirRlxStripmapReader

# Глобальные настройки стиля графиков
plt.rcParams.update({
    'font.size': 18,  # Размер шрифта 18pt
    'axes.titlesize': 18,  # Размер заголовков 18pt
    'axes.labelsize': 18,  # Размер подписей осей 18pt
    'xtick.labelsize': 16,  # Размер меток оси X
    'ytick.labelsize': 16,  # Размер меток оси Y
})

if __name__ == "__main__":
    # папка с сырыми голограммами
    # source_folder = r"E:\radar\data\хранилище Радиолокация-X\22_04_2025\raw_golograms"
    source_folder = r"E:\radar\data\хранилище Радиолокация-X\23_04_2025\raw_golograms"
    # папка, в которую будут записаны результаты анализа
    output_folder = "results_23_04_2025"

    # ищем все файлы с голограммами
    gologram_files = list()
    for root, dirs, files in os.walk(source_folder):
        for file in files:
            if file.endswith(".bin"):
                gologram_files.append(file)

    # выбираем только те файлы, для которых отсутствуют результаты обработки
    unprocessed_files = list()
    for file in gologram_files:
        file_name = os.path.basename(file)
        name = os.path.splitext(file_name)[0]
        if not os.path.exists(os.path.join(output_folder, name)):
            unprocessed_files.append(file)

    # обрабатываем
    for i, file in enumerate(unprocessed_files):
        print(f"*********** {i + 1} / {len(unprocessed_files)} ***********")
        file_name = os.path.basename(file)
        name = os.path.splitext(file_name)[0]

        # !!! все результаты сохраняем в result_folder
        result_folder = os.path.join(output_folder, name)
        os.makedirs(result_folder)

        try:
            # загружаем голограмму
            reader = NirRlxStripmapReader(dataset=name, folder=source_folder)
            gologram = reader.read_gologram(
                start_line=0, end_line=-1,
                fill_missing_lines_with_zeros=False, remove_mean=False
            )

            # Анализ голограммы

            # Задаем динамический диапазон АЦП
            dynamic_range = [-2 ** 15, 2 ** 15 - 1]

            # Задаем частоту дискретизации (720 МГц)
            Fs = 720e6

            real = np.real(gologram)
            imag = np.imag(gologram)

            def find_zero_rows(matrix):
                zero_rows = [i for i, row in enumerate(matrix) if all(x == 0 for x in row)]
                return zero_rows

            result_zeros_real = find_zero_rows(real)
            result_zeros_imag = find_zero_rows(imag)

            # Извлечение имени голограммы
            base_name = os.path.splitext(os.path.basename(file_name))[0]

            # Задаем папку для сохранения результатов
            output_dir = result_folder

            # Создаем папку, если она не существует
            os.makedirs(output_dir, exist_ok=True)

            # 1. Поиск отсчетов за пределами динамического диапазона
            rows = real.shape[0]
            rowsi = imag.shape[0]

            min_val = dynamic_range[0]
            max_val = dynamic_range[1]

            # Находим индексы выбросов
            outlier_indices = np.where((real <= min_val) | (real >= max_val))
            outlier_rows = outlier_indices[0]
            outlier_cols = outlier_indices[1]

            outlier_indicesi = np.where((imag <= min_val) | (imag >= max_val))
            outlier_rowsi = outlier_indicesi[0]
            outlier_colsi = outlier_indicesi[1]

            # Открываем файл для записи real
            output_file_real = os.path.join(output_dir, f'real_stat_{base_name}.txt')
            with open(output_file_real, 'w', encoding='utf-8') as fid:
                # Записываем заголовок
                fid.write('Анализ real голограммы\n\n')
                fid.write(f'1. Отсчеты за пределами динамического диапазона [{min_val:.2f}, {max_val:.2f}]:\n')
                fid.write('Строка\tОтсчет\tЗначение\n')

                # Записываем выбросы
                for i in range(len(outlier_rows)):
                    row = outlier_rows[i]
                    col = outlier_cols[i]
                    val = real[outlier_rows[i], outlier_cols[i]]
                    fid.write(f'{row}\t{col}\t{val:.4f}\n')
                fid.write('\n')

                # 2. Расчет количества элементов, первышающих динамический диапазон АЦП
                number_of_counts = len(outlier_rows)
                fid.write(f'2. Количество отсчетов за пределами динамического диапазона: {number_of_counts}\n')

                # 3. Расчет МО и СКО для каждой строки
                means = np.mean(real, axis=1)
                variances = np.std(real, axis=1, ddof=1)  # ddof=1 для несмещенной оценки

                # 4. Расчет МО и СКО для каждого столбца
                means_col = np.mean(real, axis=0)
                variances_col = np.std(real, axis=0, ddof=1)  # ddof=1 для несмещенной оценки

                # 5. Расчет СКО-МО и СКО-СКО real
                STD_mo_real_row = np.std(means)
                STD_sko_real_row = np.std(variances)
                STD_mo_real_col = np.std(means_col)
                STD_sko_real_col = np.std(variances_col)

                fid.write('\n3.1 CКО-МО и СКО-СКО:\n')
                fid.write('\n СКО-МО real row:\n')
                fid.write(f'{STD_mo_real_row}')
                fid.write('\n СКО-СКО real row:\n')
                fid.write(f'{STD_sko_real_row}')
                fid.write('\n СКО-МО real col:\n')
                fid.write(f'{STD_mo_real_col}')
                fid.write('\n СКО-СКО real col:\n')
                fid.write(f'{STD_sko_real_col}')

                fid.write('\n3. Математическое ожидание и СКО по строкам:\n')
                fid.write('Строка\tМат.ожидание\t\tСКО\n')
                for i in range(rows):
                    fid.write(f'{i}\t{means[i]:.4f}\t{variances[i]:.4f}\n')

                fid.write('\n4. Математическое ожидание и СКО по столбцам:\n')
                fid.write('Столбец\tМат.ожидание\t\tСКО\n')
                for i in range(real.shape[1]):
                    fid.write(f'{i}\t{means_col[i]:.4f}\t{variances_col[i]:.4f}\n')

                fid.write('\n5. Данные о МО для global_analysis:\n')
                fid.write('\n МО real-голограммы:\n')
                fid.write(f'{np.mean(real)}')
                fid.write('\n МО imag-голограммы:\n')
                fid.write(f'{np.mean(imag)}')
                fid.write('\n max МО real-голограммы (по строкам):')
                fid.write(f'{np.max(means)}')
                fid.write('\n min МО real-голограммы (по строкам):')
                fid.write(f'{np.min(means)}')
                fid.write('\n mean МО real-голограммы (по строкам):')
                fid.write(f'{np.mean(means)}')
                fid.write('\n max МО real-голограммы (по столбцам):')
                fid.write(f'{np.max(means_col)}')
                fid.write('\n min МО real-голограммы (по столбцам):')
                fid.write(f'{np.min(means_col)}')
                fid.write('\n mean МО real-голограммы (по столбцам):')
                fid.write(f'{np.mean(means_col)}')

                fid.write('\n6. Данные о СКО для global_analysis:\n')
                fid.write('\n max СКО real-голограммы (по строкам):')
                fid.write(f'{np.max(variances)}')
                fid.write('\n min СКО real-голограммы (по строкам):')
                fid.write(f'{np.min(variances)}')
                fid.write('\n mean СКО real-голограммы (по строкам):')
                fid.write(f'{np.mean(variances)}')
                fid.write('\n max СКО real-голограммы (по столбцам):')
                fid.write(f'{np.max(variances_col)}')
                fid.write('\n min СКО real-голограммы (по столбцам):')
                fid.write(f'{np.min(variances_col)}')
                fid.write('\n mean СКО real-голограммы (по столбцам):')
                fid.write(f'{np.mean(variances_col)}')

            # Первый график - Математическое ожидание (real) построчно
            plt.figure(figsize=(10, 8))
            plt.plot(range(0, rows), means)
            plt.fill_between(range(0, rows),
                             means - variances,
                             means + variances,
                             alpha=0.3, label='±СКО')
            plt.title('Математическое ожидание real по строкам')
            plt.xlabel('Номер строки')
            plt.ylabel('Мат.ожидание')
            plt.grid(True)
            plt.legend()

            # Сохраняем график
            output_img_real_row = os.path.join(output_dir, f'real_stat_row_{base_name}.png')
            plt.savefig(output_img_real_row)
            plt.close()

            # Второй график - Математическое ожидание (real) по столбцам
            plt.figure(figsize=(10, 8))
            plt.plot(range(0, real.shape[1]), means_col)
            plt.fill_between(range(0, real.shape[1]),
                             means_col - variances_col,
                             means_col + variances_col,
                             alpha=0.3, label = '±СКО')
            plt.title('Математическое ожидание real по столбцам')
            plt.xlabel('Номер столбца')
            plt.ylabel('Мат.ожидание')
            plt.grid(True)
            plt.legend()

           # Сохраняем график
            output_img_real_col = os.path.join(output_dir, f'real_stat_col_{base_name}.png')
            plt.savefig(output_img_real_col)
            plt.close()

            print(f'Анализ завершен. Результаты real сохранены в файл: {output_file_real}')
            print(f'Графики real для строк сохранены в файл: {output_img_real_row}')
            print(f'Графики real для столбцов сохранены в файл: {output_img_real_col}')

            # Открываем файл для записи imag
            output_file_imag = os.path.join(output_dir, f'imag_stat_{base_name}.txt')
            with open(output_file_imag, 'w', encoding='utf-8') as fidi:
                # Записываем заголовок
                fidi.write('Анализ imag голограммы\n\n')
                fidi.write(f'1. Отсчеты за пределами динамического диапазона [{min_val:.2f}, {max_val:.2f}]:\n')
                fidi.write('Строка\tОтсчет\tЗначение\n')

                # Записываем выбросы
                for i in range(len(outlier_rowsi)):
                    rowi = outlier_rowsi[i]
                    coli = outlier_colsi[i]
                    vali = imag[outlier_rowsi[i], outlier_colsi[i]]
                    fidi.write(f'{row}\t{col}\t{val:.4f}\n')
                fidi.write('\n')

                # 2. Расчет количества элементов, первышающих динамический диапазон АЦП
                number_of_countsi = len(outlier_rowsi);
                fidi.write(f'2. Количество отсчетов за пределами динамического диапазона: {number_of_countsi}\n')

                # 3. Расчет МО и СКО для каждой строки
                meansi = np.mean(imag, axis=1)
                variancesi = np.std(imag, axis=1, ddof=1)  # ddof=1 для несмещенной оценки

                fidi.write('\n3. Математическое ожидание и СКО по строкам:\n')
                fidi.write('Строка\tМат.ожидание\t\tСКО\n')
                for i in range(rowsi):
                    fidi.write(f'{i}\t{meansi[i]:.4f}\t{variancesi[i]:.4f}\n')

                # 4. Расчет МО и СКО для каждого столбца
                meansi_col = np.mean(imag, axis=0)
                variancesi_col = np.std(imag, axis=0, ddof=1)  # ddof=1 для несмещенной оценки

                fidi.write('\n4. Математическое ожидание и СКО по столбцам:\n')
                fidi.write('Столбец\tМат.ожидание\t\tСКО\n')
                for i in range(imag.shape[1]):
                    fidi.write(f'{i}\t{meansi_col[i]:.4f}\t{variancesi_col[i]:.4f}\n')

                # 5. Расчет СКО-МО и СКО-СКО real
                STD_mo_imag_row = np.std(meansi)
                STD_sko_imag_row = np.std(variancesi)
                STD_mo_imag_col = np.std(meansi_col)
                STD_sko_imag_col = np.std(variancesi_col)

                fidi.write('\n3.1 CКО-МО и СКО-СКО:\n')
                fidi.write('\n СКО-МО imag row:\n')
                fidi.write(f'{STD_mo_imag_row}')
                fidi.write('\n СКО-СКО imag row:\n')
                fidi.write(f'{STD_sko_imag_row}')
                fidi.write('\n СКО-МО imag col:\n')
                fidi.write(f'{STD_mo_imag_col}')
                fidi.write('\n СКО-СКО imag col:\n')
                fidi.write(f'{STD_sko_imag_col}')

                fidi.write('\n5. Данные о МО для global_analysis:\n')
                fidi.write('\n max МО imag-голограммы (по строкам):')
                fidi.write(f'{np.max(meansi)}')
                fidi.write('\n min МО imag-голограммы (по строкам):')
                fidi.write(f'{np.min(meansi)}')
                fidi.write('\n mean МО imag-голограммы (по строкам):')
                fidi.write(f'{np.mean(meansi)}')
                fidi.write('\n max МО imag-голограммы (по столбцам):')
                fidi.write(f'{np.max(meansi_col)}')
                fidi.write('\n min МО imag-голограммы (по столбцам):')
                fidi.write(f'{np.min(meansi_col)}')
                fidi.write('\n mean МО imag-голограммы (по столбцам):')
                fidi.write(f'{np.mean(meansi_col)}')

                fidi.write('\n6. Данные о СКО для global_analysis:\n')
                fidi.write('\n max СКО imag-голограммы (по строкам):')
                fidi.write(f'{np.max(variancesi)}')
                fidi.write('\n min СКО imag-голограммы (по строкам):')
                fidi.write(f'{np.min(variancesi)}')
                fidi.write('\n mean СКО imag-голограммы (по строкам):')
                fidi.write(f'{np.mean(variancesi)}')
                fidi.write('\n max СКО imag-голограммы (по столбцам):')
                fidi.write(f'{np.max(variancesi_col)}')
                fidi.write('\n min СКО imag-голограммы (по столбцам):')
                fidi.write(f'{np.min(variancesi_col)}')
                fidi.write('\n mean СКО imag-голограммы (по столбцам):')
                fidi.write(f'{np.mean(variancesi_col)}')

            # Третий график - Математическое ожидание (imag) построчно
            plt.figure(figsize=(10, 8))
            plt.plot(range(0, rowsi), meansi, color='#cc3333')
            plt.fill_between(range(0, rowsi),
                             meansi - variancesi,
                             meansi + variancesi,
                             alpha=0.3, label='±СКО', color='#cc3333')
            plt.title('Математическое ожидание imag по строкам')
            plt.xlabel('Номер строки')
            plt.ylabel('Мат.ожидание')
            plt.grid(True)
            plt.legend()

            # Сохраняем графики
            output_img_imag_row = os.path.join(output_dir, f'imag_stat_row_{base_name}.png')
            plt.savefig(output_img_imag_row)
            plt.close()

            # Четвертый график - Математическое ожидание (imag) по столбцам
            plt.figure(figsize=(10, 8))
            plt.plot(range(0, imag.shape[1]), meansi_col, color='#cc3333')
            plt.fill_between(range(0, imag.shape[1]),
                             meansi_col - variancesi_col,
                             meansi_col + variancesi_col,
                             alpha=0.3, label='±СКО', color='#cc3333')
            plt.title('Математическое ожидание imag по столбцам')
            plt.xlabel('Номер столбца')
            plt.ylabel('Мат.ожидание')
            plt.grid(True)
            plt.legend()

            # Сохраняем графики
            output_img_imag_col = os.path.join(output_dir, f'imag_stat_col_{base_name}.png')
            plt.savefig(output_img_imag_col)
            plt.close()

            print(f'Анализ завершен. Результаты imag сохранены в файл: {output_file_imag}')
            print(f'Графики imag для строк сохранены в файл: {output_img_imag_row}')
            print(f'Графики imag для столбцов сохранены в файл: {output_img_imag_col}')

            ## Спектральный анализ по дальности для комплексного сигнала

            # Инициализация массива для накопления спектра

            for i in range(rows - 1):
                # 1. Формирование комплексного сигнала
                complex_line = real[i, :] + 1j * imag[i, :]

                # 2. Расчет спектра комплексного сигнала
                spectrum = np.fft.fftshift(np.abs(np.fft.fft(complex_line)))
                spectrum_db = 10 * np.log10(spectrum)  # Добавляем малое значение для избежания log(0)

                # 3. Частотная ось
                N = len(complex_line)  # количество отсчетов в строке
                freqs = np.fft.fftfreq(N, d=1 / Fs)  # частотная ось в Гц
                freqs_shifted = np.fft.fftshift(freqs) / 1e6  # сдвиг нуля и перевод в МГц

                # 4. Суммирование
                spectrum_db += spectrum_db

            # 5. Усреднение
            spectrum_avg = spectrum_db / rows

            # Визуализация усредненного спектра по дальности
            plt.figure(figsize=(10, 8))
            plt.plot(freqs_shifted, spectrum_avg)
            plt.title('Усредненный спектр сигнала по дальности')
            plt.xlabel('Частота, МГц')
            plt.ylabel('Мощность, дБ')
            plt.grid(True)

            # 10. Сохранение графика
            output_spectrum = os.path.join(output_dir, f'spectrum_{base_name}.png')
            plt.tight_layout()
            plt.savefig(output_spectrum)
            plt.close()

            print(f"Спектральный анализ по дальности завершен. График сохранен: {output_spectrum}")

            ## Расчет мат.ожидания для real и мнимой составляющей (вся голограмма)
            mean_real_gologram = np.mean(real)
            mean_imag_gologram = np.mean(imag)

            # Содание Typst-файла
            typ_filename = os.path.join(output_dir, f"analysis_{base_name}.typ")

            # Экранируем специальные символы в имени файла
            safe_base_name = base_name.replace("{", "\\{").replace("}", "\\}").replace("_", "\\_")

            typ_content = f"""
            #set page(width: auto, height: auto, margin: 1cm)
            #align(center)[#text(size: 28pt, weight: "bold")[Анализ голограммы]]
            #align(center)[#text(size: 28pt)[{safe_base_name}]
            ]
    
            #figure(
        table(
            columns: 2,
            align: center,
            stroke: (x: 0.5pt, y: 0.5pt),
        
            // Заголовки столбцов
            [#text(size: 22pt, weight: "bold")[Математическое ожидание голограммы]], 
            [#text(size: 22pt, weight: "bold")[Значение]],
            
            // Строки данных
            [#text(size: 22pt)[МО real-голограммы]], 
            [#text(size: 22pt)[{int(np.round(mean_real_gologram))}]],
            
            [#text(size: 22pt)[МО imag-голограммы]], 
            [#text(size: 22pt)[{int(np.round(mean_imag_gologram))}]],
                           
            [#text(size: 22pt)[Количество отсчетов за пределами ДД (real)]], 
            [#text(size: 22pt)[{int(number_of_counts)}]],
                        
            [#text(size: 22pt)[Количество отсчетов за пределами ДД (imag)]], 
            [#text(size: 22pt)[{int(number_of_countsi)}]],
                        
            [#text(size: 22pt)[Номера нулевых строк в real]], 
            [#text(size: 22pt)[{result_zeros_real}]],
                        
            [#text(size: 22pt)[Номера нулевых строк в imag]], 
            [#text(size: 22pt)[{result_zeros_imag}]]
                ),
            )
    
            #align(center)[ #text(size: 26pt, weight: "bold")[Анализ мат.ожидания и СКО]]
            
            #figure(
                grid(
                    columns: 2,
                    gutter: 2mm,
                    [#image("real_stat_row_{base_name}.png")], 
                    [#image("real_stat_col_{base_name}.png")],
                ),
            )
 
            #figure(
                grid(
                    columns: 2,
                    gutter: 2mm,
                    [#image("imag_stat_row_{base_name}.png")], 
                    [#image("imag_stat_col_{base_name}.png")],
                ),
            )
    
            #align(center)[#text(size: 26pt, weight: "bold")[Спектральный анализ по дальности]]
            #align(center)[#image("spectrum_{base_name}.png")]
            """

            with open(typ_filename, 'w', encoding='utf-8') as typ_file:
                typ_file.write(typ_content)

            # Сохраняем typ-файл
            typ_path = os.path.join(output_dir, f"analysis_{base_name}.typ")
            pdf_path = typ_path.replace('.typ', '.pdf')

            # Вызываем компилятор Typst
            typst_path = r'C:\Users\Paneeva\Documents\libraries\typst\typst.exe'

            result = subprocess.run([typst_path, 'compile', typ_path, pdf_path])

            # Удаляем все лишние файлы в папке, чтобы память не занимали
            files_del = [f"real_stat_row_{base_name}.png",
                         f"real_stat_col_{base_name}.png",
                         f"imag_stat_row_{base_name}.png",
                         f"imag_stat_col_{base_name}.png",
                         f"spectrum_{base_name}.png",
                         f"analysis_{base_name}.typ"]
            for files_name in files_del:
                file_path = os.path.join(result_folder, files_name)
                if os.path.exists(file_path):
                    os.remove(file_path)
            print(f"Все временные файлы удалены")
        except BaseException:
            shutil.rmtree(result_folder)
            raise






