import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from cft import cft, variance_spectrum

st.set_page_config(page_title="Fourier Clásica", page_icon="〰️", layout="centered")

st.title("〰️ Transformada Clásica de Fourier")
st.markdown(
    "Arma una señal sumando armónicos y ruido, y mira su descomposición "
    "en frecuencias con la **Transformada Clásica de Fourier** "
    "([Bendat & Piersol, 2011](https://www.wiley.com/en-us/Random+Data%3A+Analysis+and+Measurement+Procedures%2C+4th+Edition-p-9780470248775)), "
    "el mismo método usado en el curso de Series de Tiempo y Análisis "
    "Espectral."
)

with st.sidebar:
    st.header("Señal")
    N = st.slider("Largo de la serie (N)", 100, 2000, 1000, 50)
    dt = 1.0
    n_arm = st.slider("Número de armónicos", 1, 5, 3)

    periods, amps = [], []
    for i in range(n_arm):
        c1, c2 = st.columns(2)
        periods.append(c1.slider(f"Periodo {i+1}", 4, 200, 10 * (i + 1) ** 2, key=f"p{i}"))
        amps.append(c2.slider(f"Amplitud {i+1}", 0.1, 2.0, 0.5, 0.1, key=f"a{i}"))

    ruido = st.slider("Desviación estándar del ruido", 0.0, 3.0, 0.5, 0.1)
    semilla = st.number_input("Semilla aleatoria", 0, 9999, 0)

    st.header("Visualización")
    n_mostrar = st.slider("Armónicos individuales a mostrar", 0, n_arm, min(n_arm, 3))

t = np.arange(1, N + 1) * dt
armonicos = [amps[i] * np.sin(2 * np.pi * t / periods[i]) for i in range(n_arm)]
S = np.sum(armonicos, axis=0)
rng = np.random.default_rng(semilla)
x = S + rng.normal(0, ruido, size=N) if ruido > 0 else S.copy()

A0, Aq, Bq = cft(x)
ev = variance_spectrum(Aq, Bq)
freqs = np.arange(1, len(ev) + 1) / (N * dt)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6))

ax1.plot(t, x, color="#2c3e50", lw=0.9)
ax1.set_xlabel("Tiempo")
ax1.set_ylabel("Amplitud")
ax1.set_title("Serie sintética")
ax1.grid(alpha=0.3)

ax2.plot(freqs, ev, color="#c0392b", lw=1)
for p in periods:
    ax2.axvline(1 / p, color="0.6", lw=1, ls="--")
ax2.set_xlabel("Frecuencia [1/dt]")
ax2.set_ylabel(r"Varianza [$c_q^2/2$]")
ax2.set_title("Espectro de varianza (líneas punteadas = frecuencias reales)")
ax2.set_xlim(0, max(2 / min(periods), 0.05))
ax2.grid(alpha=0.3)

fig.tight_layout()
st.pyplot(fig)
plt.close(fig)

col1, col2, col3 = st.columns(3)
col1.metric("Varianza de x", f"{np.var(x):.3f}")
col2.metric("Σ espectro (Parseval)", f"{ev.sum():.3f}")
col3.metric("Armónico dominante", f"{freqs[np.argmax(ev)]:.4f} [1/dt]")

if n_mostrar > 0:
    st.subheader(f"Primeros {n_mostrar} armónicos que componen la serie")
    fig3, axes = plt.subplots(n_mostrar, 1, figsize=(7, 1.8 * n_mostrar), sharex=True)
    axes = np.atleast_1d(axes)
    for i in range(n_mostrar):
        axes[i].plot(t, armonicos[i], color="#2980b9", lw=1)
        axes[i].set_ylabel("Amplitud")
        axes[i].set_title(f"Armónico {i+1} — periodo {periods[i]}, amplitud {amps[i]}", fontsize=10)
        axes[i].grid(alpha=0.3)
    axes[-1].set_xlabel("Tiempo")
    fig3.tight_layout()
    st.pyplot(fig3)
    plt.close(fig3)

with st.expander("¿Cómo se calcula?"):
    st.markdown(
        r"""
Para una serie $x_n$ ($n=1,\dots,N$) con paso $\Delta t$ constante:

$$x_n = A_0 + \sum_{q=1}^{Q}\left[A_q\cos\!\left(\frac{2\pi q n}{N}\right) + B_q\sin\!\left(\frac{2\pi q n}{N}\right)\right]$$

con $A_q=\frac{2}{N}\sum_n x_n\cos(\cdot)$, $B_q=\frac{2}{N}\sum_n x_n\sin(\cdot)$,
y el espectro de varianza por armónico $(A_q^2+B_q^2)/2$, que por el
**teorema de Parseval** suma exactamente la varianza total de la serie
(confirmado arriba: "Varianza de x" ≈ "Σ espectro").

Esta es una implementación en Python del `cft.m` original del curso, con
una corrección: cuando $N$ es par existe un término extra de Nyquist
($q=N/2$), que en el `.m` original sobrescribe el último elemento en vez
de agregarse, y al que luego se le aplicaba por error el mismo factor
$1/2$ del resto del espectro — acá se agrega como componente aparte y se
contabiliza sin ese factor, tal como exige Parseval. Verificado con
reconstrucción exacta, Parseval, y contraste contra la FFT de NumPy en
[`tests/test_cft.py`](https://github.com/IPartarrieu/Classic-Fourier-Transformation/blob/main/fourier/tests/test_cft.py).
"""
    )
