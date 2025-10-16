import numpy as np

def generate_sinc_signal(time_range, f_discr, sinc_width):
    t = np.linspace(-time_range, time_range, f_discr)
    sinc_linear = np.sinc(t / sinc_width)
    sinc_norm = np.abs(sinc_linear) / np.max(np.abs(sinc_linear))
    sinc_db = 20 * np.log10(sinc_norm)

    return t, sinc_db, sinc_norm