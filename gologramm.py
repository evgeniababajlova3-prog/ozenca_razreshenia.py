import h5py
import numpy as np


def load_radar_image_from_hdf5(hdf5_path):
    """
    Прямое чтение голограммы из HDF5 файла через h5py
    """
    try:
        with h5py.File(hdf5_path, 'r') as f:
            print("Структура HDF5 файла:")

            # Функция для вывода структуры файла
            def print_structure(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  Dataset: {name}, shape: {obj.shape}, dtype: {obj.dtype}")
                elif isinstance(obj, h5py.Group):
                    print(f"  Group: {name}")

            f.visititems(print_structure)

            # Пытаемся найти данные разными способами
            radar_image_complex = None

            # Способ 1: Ищем готовый комплексный массив
            for name in f:
                if isinstance(f[name], h5py.Dataset):
                    dataset = f[name]
                    if dataset.dtype == np.complex64 or dataset.dtype == np.complex128:
                        radar_image_complex = dataset[:]
                        print(f"Найден комплексный массив: {name}")
                        break

            # Способ 2: Ищем отдельно real и imag
            if radar_image_complex is None:
                real_data = None
                imag_data = None

                # Ищем real часть
                for name in f:
                    if isinstance(f[name], h5py.Dataset):
                        if 'real' in name.lower() or 're' in name.lower():
                            real_data = f[name][:]
                            print(f"Найден real: {name}")
                        elif 'imag' in name.lower() or 'im' in name.lower():
                            imag_data = f[name][:]
                            print(f"Найден imag: {name}")

                if real_data is not None and imag_data is not None:
                    radar_image_complex = real_data + 1j * imag_data
                elif real_data is not None:
                    # Если есть только real, используем его как амплитуду
                    radar_image_complex = real_data.astype(np.complex64)
                    print("Используется только real часть")
                else:
                    # Способ 3: Берем первый попавшийся 2D массив
                    for name in f:
                        if isinstance(f[name], h5py.Dataset) and len(f[name].shape) >= 2:
                            radar_image_complex = f[name][:].astype(np.complex64)
                            print(f"Используется массив: {name}")
                            break

            if radar_image_complex is None:
                raise ValueError("Не удалось найти подходящие данные в HDF5 файле")

            # Если данные 3D, берем первый слой
            if len(radar_image_complex.shape) == 3:
                radar_image_complex = radar_image_complex[0, :, :]
                print(f"Взят первый слой из 3D данных, новый shape: {radar_image_complex.shape}")

            radar_image = np.abs(radar_image_complex)

            print(f"Загружена голограмма размером: {radar_image.shape}")
            print(f"Диапазон амплитуд: {np.min(radar_image):.6f} - {np.max(radar_image):.6f}")

            return radar_image_complex, radar_image, -80

    except Exception as e:
        print(f"Ошибка при загрузке HDF5 файла: {e}")
        raise