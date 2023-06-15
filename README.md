# Classic-Fourier-Transformation
Función para calcular la transformada de Fourirer en Matlab. La función tiene como entrada una serie de tiempo estacionaria (x) con un paso de tiempo \delta-t constante (paso de tiempo de la serie x). La salida nos entrega las componentes A_0, A_q y B_q como vectores columna.
La función está basada en la Transformada Clásica de Fourier expuesta en el libro Random data (2011).

El código cft_ejemplo.m corresponde a un ejercicio del curso Series de Tiempo y Análisis Espectral de Magister en Geofísica. Consta en la creación de una serie estacional con 5 armónicos además de ruido para simular datos reales. Se calculan las varianzas mediante un análisis tras una transformada de Fourier con la función Clasic Fourier Transformation (cft.m) para cada armónico.
El código welch_metod.m también consiste en un ejercicio, donde se calcula la densidad del espectro mediante el método de Welch, aplicando ventanas de Hanning a cada armónico.
 
Bendat, J. S., \& Piersol, A. G. (2011). Random data: analysis and measurement procedures. John Wiley & Sons.
