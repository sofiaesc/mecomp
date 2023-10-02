% PoissonDifReaccion.m
% Ecuacion de Poisson con condiciones Robin
%
% Programa que resuelve la ecuacion de Poisson con condiciones de contorno
% Dirichlet-Robin, usando Diferencias Finitas
% Calcula temperatura en el extremo derecho y flujo en el izquierdo
% Calcula la ecuacion de difusión con  condiciones Dirichlet y Mixtas
%
% -K0*u''(x)+ pR*u = Q(x), 0<x<L
% u(0) = a; K0*u'(L) + H1*u(L) = H2*uE

% paramwetros del problema:
L=2; K0=0.57; a=6;
pR=@(x) 0.1*x.^3+2.5;
H1=15; H2=15; uE=4;

% FUENTE DEL PROBLEMA
%Q =@(x) (100*exp(-20*(x-0.3).^2));
%Q=@(x) (100*(x<0.7));
%Q=@(x) 20*sin(5*(x-1));
Q = @(x) 2*x.*(2-x);


% parametros del metodo de resolucion
N=50; % SUBINTERVALOES
x = linspace(0,L,N+2)';
size(x)

%% Armado de la matriz
h = L/N;
unos = ones(N+2,1);
size(unos)
diagonales = [-1*unos (2+((h.^2/K0).*pR(x(1:N+2)))).*unos -1*unos];
matriz = spdiags(diagonales, [-1 0 1], N+2, N+2);
matriz(1,[1:2])=[1 0];
matriz(N+2, [N:N+2]) = [-1 2*h*H1/K0 1];

% Armado del lado derecho
F =[a; (h^2/K0)*Q(x(2:N+1)); 2*h*H2*uE/K0];
size(F)

% Resolucion
U = matriz\F;
size(U)
U(N+2) = []; % otra forma: U = U(1:N+1);
% Otra forma
%
xx=x(1:N+1);
figure(1)
plot(xx,U,'*-')
grid on;
grid minor;

% temperatura en el extremo derecho
TexD=U(N+1);
disp('Temperatura en el extremo derecho:');
disp(TexD)

FlujoExtIzq = (U(1)-U(2))/h;
disp('Flujo Extremo Izquierdo:');
disp(FlujoExtIzq)

FlujoExtDer=H1*(U(N+1)-uE)
