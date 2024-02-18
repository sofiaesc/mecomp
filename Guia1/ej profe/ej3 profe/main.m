close all; clear all;

 # Malla:
 L1=0; L2=1;
 N = 40;
 dx = (L2-L1)/N;
 xnode = L1:dx:L2;

 # Condiciones de borde:
 cb = [1 0 -1;2 0 -1];

 # Constantes del modelo:
 model.rho = 1;
 model.cp = 1;
 model.k = 1;
 model.c = 0;

 # Parte no estacionaria:
 lambda = 4;
 alpha = model.k/(model.rho*model.cp);
 dt = lambda*(dx*dx)/alpha;
 model.G = 100*xnode.^0;
 et = [2 1/dt 1e-12 dt];

 # T anterior
 Tant = -10*(xnode.*xnode)+20*xnode;

 # Solución aproximada:
 T = difFinitas(xnode, model, cb, et, Tant);

hold on
plot(xnode, Tex, 'r')
xlabel('x'); ylabel('Temp')
legend('Aproximada', 'Exacta')
