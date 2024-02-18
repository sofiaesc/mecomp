# ESTE SCRIPT ES ESCRIBIENDO LA MATRIZ QUE OBTUVIMOS EN EL DESARROLLO A PATA
clf; clear all; clc;

# DATOS DE LA FUNCIÓN
phi = @(x) exp(-0.1*x.^2);          # Función
dphi = @(x) -0.2*x*exp(-0.1*x.^2);       # Derivada de la función
i = -0.2;                  # Punto
real = dphi(i);           # Solución analítica

puntos = [-2 -1 0];      # Qué puntos quiero usar
# Pongo el sistema de ecuaciones que encontré acá.
K = [1 1 1;
     -7/4 -3/4 0;
     49/16 9/16 0];

# Pasos de tiempo:
n = 10;
h = 1./(2.^(0:n-1));
p = [1 2];                    # Puntos intermedios a tomar para calcular p
errores = [];

for j=p(1):p(2)
  f = [0; 1/h(j); 0];            # Defino la F que encontré
  coef = K\f;
  coords = i + puntos*h(j);          # Obtengo coordenada de cada punto con esa h
  val = phi(coords);                 # Evalúo el punto con la función phi
  aprox = dot(coef,val);             # Calculo la aproximación
  errores = [errores, abs(real - aprox)];
endfor

h = h(p(1):p(2));
plot(log(h),log(errores))
poli = polyfit(log(h),log(errores),1);
p = poli(1)                          # orden
