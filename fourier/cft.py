"""Transformada Clasica de Fourier (CFT) para series de tiempo estacionarias.

Puerto a Python del `cft.m` original (curso Series de Tiempo y Analisis
Espectral, Magister en Geofisica), basado en Bendat & Piersol (2011)
"Random Data": para una serie x de largo N con paso de tiempo dt
constante,

    x_n = A0 + sum_{q=1}^{Q} [Aq*cos(2*pi*q*n/N) + Bq*sin(2*pi*q*n/N)]

Diferencia respecto al .m original: en `cft.m`, el termino de Nyquist
(q = N/2, solo cuando N es par) se calcula con la formula correcta
(normalizacion 1/N en vez de 2/N, sin componente seno) pero
SOBRESCRIBE el ultimo elemento de Aq en lugar de agregarse como un
componente nuevo, y luego el script de ejemplo (cft_ejemplo.m) le
aplica la misma formula de varianza (Aq^2+Bq^2)/2 que a los demas
terminos, cuando por el teorema de Parseval el termino de Nyquist no
lleva el factor 1/2 (su base cos(n*pi) tiene norma llena, no la mitad,
porque no oscila en amplitud). Aca se implementa la version correcta:
el Nyquist se agrega como componente adicional y se contabiliza aparte
en el espectro de varianza. Ver tests/test_cft.py para la verificacion
(identidad de reconstruccion exacta y teorema de Parseval).
"""

import numpy as np


def cft(x):
    """Calcula A0, Aq, Bq de la serie x (array 1D).

    Returns
    -------
    A0 : float
    Aq : ndarray, forma (Q,) donde Q = N//2 (incluye el termino de
        Nyquist como ultimo elemento si N es par).
    Bq : ndarray, forma (Q,) si N es impar, o (Q-1,) si N es par (el
        termino de Nyquist no tiene componente seno).
    """
    x = np.asarray(x, dtype=float)
    N = len(x)
    n = np.arange(1, N + 1)
    A0 = float(x.mean())

    has_nyquist = N % 2 == 0
    Q = N // 2
    Q_ac = Q - 1 if has_nyquist else Q

    q = np.arange(1, Q_ac + 1).reshape(-1, 1)
    theta = (2 * np.pi / N) * q * n.reshape(1, -1)
    Aq = (2 / N) * (np.cos(theta) @ x)
    Bq = (2 / N) * (np.sin(theta) @ x)

    if has_nyquist:
        a_nyquist = (np.cos(n * np.pi) @ x) / N
        Aq = np.concatenate([Aq, [a_nyquist]])

    return A0, Aq, Bq


def variance_spectrum(Aq, Bq):
    """Espectro de varianza (Aq^2+Bq^2)/2 por armonico general, con el
    termino de Nyquist (si existe, cuando Aq tiene un elemento mas que
    Bq) contabilizado sin el factor 1/2, segun el teorema de Parseval.
    """
    Q_ac = len(Bq)
    ev = (Aq[:Q_ac] ** 2 + Bq**2) / 2
    if len(Aq) > Q_ac:
        ev = np.concatenate([ev, [Aq[-1] ** 2]])
    return ev


def reconstruct(A0, Aq, Bq, N):
    """Reconstruye la serie original de largo N a partir de A0, Aq, Bq."""
    n = np.arange(1, N + 1)
    Q_ac = len(Bq)
    q = np.arange(1, Q_ac + 1).reshape(-1, 1)
    theta = (2 * np.pi / N) * q * n.reshape(1, -1)
    x = A0 + Aq[:Q_ac] @ np.cos(theta) + Bq @ np.sin(theta)
    if len(Aq) > Q_ac:
        x = x + Aq[-1] * np.cos(n * np.pi)
    return x
