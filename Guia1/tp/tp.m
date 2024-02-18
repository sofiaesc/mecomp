clc; close all; clear all;

 #--- Parámetros para la función de diferencias finitas ---#

 L1=0; L2=1;                                #Extremos laterales
 tini=0; tfin=4;                            #Intervalo de tiempo
 N = 40;                                    #Cantidad de intervalos
 dx = (L2-L1)/N;                            #Delta x
 xnode = L1:dx:L2;                          #Vector de nodos
 cb = [1 0 -1; 2 0 -1];                     #Condiciones en los extremos
 model.rho = 1; model.cp = 1;               #Rho y Cp
 model.k = 1; model.c = 0;                  #K y C
 model.G = 100.*xnode.^0;                   #Fuente hasta t=2

 dx = max(diff(xnode));
 alpha = model.k/(model.rho*model.cp);
 #dt = 0.9*(0.5*(dx*dx))/alpha;             #dt para el esquema explicito
 dt = 4*(dx*dx)/alpha;                      #dt para el esquema implicito

 it = round(4/dt);                          #Numero de iteraciones segun dt
 tol = 1e-100;                              #Tolerancia para el corte
 mode = 3;                                  # estacionario, 1: explicito, 2: implicito, 3: semi (CN)
 et = [mode it tol dt];

 Tinicial = (-10.*xnode.^2 + 20.*xnode)';   #Temperatura en t=0

 T = tp_dif_finitas(xnode, model, cb, et, Tinicial);
 plot(xnode, Tinicial, 'r', 'linewidth', 1.5) #Gráfica para Tinicial
 temp_max = max(max(T));
 axis([0, 1, 0, temp_max]);
 xlabel('Posición en la barra'); ylabel('Temperatura'); hold on;

 #--- Animación evolución Temperatura en todos los nodos ---#

 for i=1:round(it/100):it
     if(i*dt>=2)
          plot(xnode, T(:,i), '-b')
          title(sprintf('nit: %d - tiempo: %f',i,i*dt));
          hold on
          pause(0.05)
        else
          plot(xnode, T(:,i), '-g')
          title(sprintf('nit: %d - tiempo: %f',i,i*dt));
          hold on
          pause(0.05)
        endif
  endfor

  nodo = N/2 + 1;

 #--- Gráfica evolución de Temperatura en el punto medio ---#
    figure(2)
    xlim([tini,tfin]);
    hold on
    x_aux = 0:et(4):(it*et(4));
    plot(x_aux,T(nodo,:),'r', 'linewidth', 2); hold on;
    plot(x_aux,7.*(x_aux.^0),'b', 'linewidth', 2); hold on;
    title("Temperatura en el punto medio");
    xlabel('Tiempo'); ylabel('Temperatura');

#--- Instante de tiempo donde la temperatura baja de 7° en el punto medio ---#
    for i=1:it
      if(T(nodo,i) >= 7 && T(nodo,i+1) < 7)
        printf("Temperatura en el punto medio para el instante %f: %f \n", i*et(4), T(nodo,i))
        printf("Temperatura en el punto medio para el instante %f: %f \n", (i+1)*et(4), T(nodo,i+1))
      endif
    endfor
