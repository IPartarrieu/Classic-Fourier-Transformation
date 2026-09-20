import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cft import cft, reconstruct, variance_spectrum  # noqa: E402

rng = np.random.default_rng(0)

# 1) Coseno puro: toda la energia debe caer en el bin q=k, con Aq[k-1]=C
for N in [100, 101]:  # par e impar
    k, C = 7, 2.5
    n = np.arange(1, N + 1)
    x = C * np.cos(2 * np.pi * k * n / N)
    A0, Aq, Bq = cft(x)
    assert abs(A0) < 1e-8
    assert abs(Aq[k - 1] - C) < 1e-8, (N, Aq[k - 1])
    assert abs(Bq[k - 1]) < 1e-8
    other_Aq = np.delete(Aq[: len(Bq)], k - 1)
    other_Bq = np.delete(Bq, k - 1)
    assert np.allclose(other_Aq, 0, atol=1e-8)
    assert np.allclose(other_Bq, 0, atol=1e-8)

# 2) Seno puro: analogo, energia en Bq[k-1]
N, k, C = 200, 15, 1.7
n = np.arange(1, N + 1)
x = C * np.sin(2 * np.pi * k * n / N)
A0, Aq, Bq = cft(x)
assert abs(Bq[k - 1] - C) < 1e-8
assert abs(Aq[k - 1]) < 1e-8

# 3) Reconstruccion exacta (N par e impar) para una serie ruidosa cualquiera
for N in [50, 51, 300, 301]:
    x = rng.normal(size=N)
    A0, Aq, Bq = cft(x)
    recon = reconstruct(A0, Aq, Bq, N)
    assert np.allclose(recon, x, atol=1e-8), N

# 4) Teorema de Parseval: varianza de x == suma del espectro de varianza
for N in [80, 81]:
    x = rng.normal(size=N) * 3 + 5
    A0, Aq, Bq = cft(x)
    ev = variance_spectrum(Aq, Bq)
    assert abs(np.var(x) - ev.sum()) < 1e-8, (N, np.var(x), ev.sum())

# 5) Contraste contra la FFT de numpy (referencia externa confiable).
#    cft() usa n=1..N (como el cft.m original), mientras que rfft usa la
#    convencion estandar n=0..N-1; eso rota la fase de cada Aq/Bq pero no
#    cambia la magnitud (energia) de cada armonico, que es la cantidad que
#    importa para el espectro de varianza. Se compara magnitud, no fase.
N = 128
x = rng.normal(size=N)
A0, Aq, Bq = cft(x)
X = np.fft.rfft(x)  # X[q] = sum_n x_n * exp(-2pi*i*q*n/N), n=0..N-1 (0-indexed)
mag_cft = Aq[: N // 2 - 1] ** 2 + Bq**2
mag_fft = ((2 / N) * X.real[1 : N // 2]) ** 2 + ((2 / N) * X.imag[1 : N // 2]) ** 2
assert np.allclose(mag_cft, mag_fft, atol=1e-8)
# Nyquist (bin N/2 del rfft) puramente real, con normalizacion 1/N:
assert np.allclose(Aq[-1] ** 2, (X.real[N // 2] / N) ** 2, atol=1e-8)
assert abs(X.imag[N // 2]) < 1e-8  # el bin de Nyquist siempre es real

print("Todas las verificaciones de la transformada clasica de Fourier pasaron OK")
