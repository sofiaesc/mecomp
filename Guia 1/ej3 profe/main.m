close all; clear all;

% L1 y L2: extremos izquierdo y derecho respectivamente
% N: cantidad de puntos
% cb: condiciones de borde [tipo valor1 valor2]
  % tipo==1 -> Dirichet -- valor1 = Temp. impuesta
  % tipo==2 -> Neumann -- valor1 = q (flujo impuesto)
  % tipo==3 -> Robin -- valor1 = h (conveccion) -- valor2 = Temp. ext.
% model: struct con los parametros del modelo incluida la fuente
%        punto a punto
% et: esquema temporal [tipo val1 val2 val3]
  % tipo==0 -> Estacionario
  % tipo==1 -> Explicito
  % tipo==2 -> Implicito
  % Si tipo==1 o tipo==2: val1=maxIt, val2=tolError
  % Si tipo==2 val3=dt

% Cualquier ejercicio puede ser resuelto en su versión transiente o estacionaria

# EJERCICIO 2.a
 L1=0; L2=1;
 N = 100;
 dx = (L2-L1)/N;
 xnode = L1:dx:L2;
 cb = [1 10 -1;1 50 -1];
 model.rho = 1; model.cp = 1;
 model.k = 1; model.c = 0;
 model.G = @(x,t) (t < 2)*100; et = [0 -1 -1 -1];
 Tex = -25*(xnode.*xnode)+65*xnode+10;
 T = difFinitas(xnode, model, cb, et);

hold on
plot(xnode, Tex, 'r')
xlabel('x'); ylabel('Temp')
legend('Aproximada', 'Exacta')
