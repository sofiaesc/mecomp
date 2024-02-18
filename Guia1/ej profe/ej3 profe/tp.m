clc; close all; clear all;

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

 L1=0; L2=1;
 N = 40;
 dx = (L2-L1)/N;
 xnode = L1:dx:L2;
 cb = [1 0 -1;2 0 -1];
 model.rho = 1; model.cp = 1;
 model.k = 1; model.c = 0;
 # Como N=40 dx = 1/40 = 0.025
 # alpha = k/rhoCp = 1/1
 dt = 4.*dx.*dx;

 model.G = 100.*xnode.^0; et = [2 2000 1e-10000 dt];
 Tant = -10.*xnode.^2 + 20.*xnode;

 axis = [0 1 0 50]
 T = tp_dif_finitas(xnode, model, cb, et,Tant');

hold on
#plot(xnode, T, 'b')
#plot(xnode, Tant, 'r')
#xlabel('x'); ylabel('Temp')
#legend('Aproximada', 'Exacta')
