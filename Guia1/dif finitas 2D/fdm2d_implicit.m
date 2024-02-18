function [PHI, Q] = fdm2d_implicit(K,F,xnode,neighb,model,dt)
% DescripciÃ³n: mÃ³dulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal implÃ­cito.

% Entrada:
% * K: matriz del sistema (difusiÃ³n + reacciÃ³n)
% * F: vector de flujo tÃ©rmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo 
% de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numÃ©rico, etc.)
% * dt: paso temporal crÃ­tico para mÃ©todo explÃ­cito.

% Salida:
% * PHI: matriz soluciÃ³n. Cada elemento del vector representa un valor escalar 
%   asociado a cada nodo de la malla, y su posiciÃ³n dentro del vector depende de
%   cÃ³mo se especificÃ³ cada nodo en xnode. Cada columna representa una iteraciÃ³n
%   del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas 
%   representa una iteraciÃ³n del esquema temporal (en total 2Ã—nit columnas).
% ----------------------------------------------------------------------
    disp('Se inicializa al calculo del modelo IMPLICITO...');
    nt = model.maxit;
    N = size(xnode,1);
    tol = model.tol;
    PHI = model.PHI_n;
    PHIn = model.PHI_n;
    Q = fdm2d_flux(PHI,neighb,xnode,model.k);
    a = (model.rho*model.cp)/dt;
    I = eye(N,N);
    Ks = a*I + K;
    for i=1:nt
        B = F + a*PHIn;
        PHIa = Ks\B;
        err = norm(PHIa-PHIn,2)/norm(PHIa,2);
        PHIn = PHIa;
        PHI = [PHI PHIa];      
        Qa = fdm2d_flux(PHIa,neighb,xnode,model.k);
        Q = [Q Qa];
        if (err < tol)
          disp('Ha completado el calculo IMPLICITO.');
          printf("\n Cantidad de iteracion: ");
          i
            return;
        endif
    endfor
    disp('Llego la maxima iteracion IMPLICITO.');
endfunction