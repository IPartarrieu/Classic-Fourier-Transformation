# Classic-Fourier-Transformation

Demo interactiva de la **Transformada Clásica de Fourier** (Bendat & Piersol, 2011, *Random Data*): arma una señal sumando armónicos y ruido, y mira su descomposición en frecuencias en tiempo real.

🔗 **Demo en vivo:** [classic-fourier-tr-nvubqfatfs2ra5ty2mqfhs.streamlit.app](https://classic-fourier-tr-nvubqfatfs2ra5ty2mqfhs.streamlit.app/)

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://classic-fourier-tr-nvubqfatfs2ra5ty2mqfhs.streamlit.app/)

## Qué hace

Elige el número de armónicos, su periodo y amplitud, y el nivel de ruido — la app recalcula al instante la serie sintética y su espectro de varianza, marcando con líneas punteadas dónde deberían caer los picos según las frecuencias reales que elegiste. También confirma en vivo el teorema de Parseval (la varianza de la señal debe ser igual a la suma del espectro).

## Cómo está construido

- **`fourier/cft.py`** — puerto a Python de la función original `cft.m` (ver [`matlab_exercise/cft.m`](matlab_exercise/cft.m)), con una corrección: el término de Nyquist (cuando N es par) se agrega como componente aparte en vez de sobrescribir el último armónico general, y se contabiliza en el espectro de varianza sin el factor 1/2 que sí llevan los demás términos, tal como exige Parseval.
- Verificado en [`fourier/tests/test_cft.py`](fourier/tests/test_cft.py): reconstrucción exacta de la señal, teorema de Parseval, y contraste de magnitudes contra la FFT de NumPy.
- **`fourier/app.py`** — interfaz interactiva en Streamlit.

## Correr en local

```bash
cd fourier
pip install -r requirements.txt
streamlit run app.py
```

## Ejercicios originales (MATLAB)

`cft.m` (la función), `cft_ejemplo.m` (armónicos + ruido) y `welch_metod.m` (densidad espectral por el método de Welch, con ventanas de Hanning) del curso Series de Tiempo y Análisis Espectral, Magíster en Geofísica, quedaron en [`matlab_exercise/`](matlab_exercise/).

Bendat, J. S., & Piersol, A. G. (2011). *Random data: analysis and measurement procedures*. John Wiley & Sons.
