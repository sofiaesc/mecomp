# ESTE SCRIPT USA LA FUNCIÓN TAYLORDF QUE ME CALCULA LOS COEFICIENTES.
clf; clear all; clc;
pkg load symbolic

phi = @(x) exp(2*x);          # Función
dphi = @(x) 2*exp(2*x);       # Derivada k de la función (CAMBIAR SEGÚN K)
i = 1.5;                  # Punto
real = dphi(i);           # Solución analítica

puntos = [-2 -1 0];      # Qué puntos quiero usar
k = 1;                    # Derivada a aproximar
[coef] = taylorDFsyms(puntos,k)

# Pasos de tiempo:
n = 10;
h = (1./(2.^(0:n-1)))/4;
p = [3 4];                    # Puntos intermedios a tomar para calcular p
errores = [];

for j=p(1):p(2)
  [coef] = taylorDF(puntos,k,h(j));  # Calculo coeficiente con función
  coords = i + puntos*h(j);          # Obtengo coordenada de cada punto con esa h
  val = phi(coords);                 # Evalúo el punto con la función phi
  aprox = dot(coef,val);             # Calculo la aproximación
  errores = [errores, abs(real - aprox)];
endfor

h = h(p(1):p(2));
plot(log(h),log(errores))
poli = polyfit(log(h),log(errores),1);
p = poli(1)                          # orden
