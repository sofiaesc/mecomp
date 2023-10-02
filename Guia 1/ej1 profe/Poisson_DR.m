clear all; clc;
% Ecuacion de Poisson con condiciones Robin
%
% Programa que resuelve la ecuacion de Poisson con condiciones de contorno
% Dirichlet-Robin, usando Diferencias Finitas
%
% -K0*u''(x) = Q(x), 0<x<L
% u(0) = a; K0*u'(L) + H1*u(L) = H2*uE


% parametros del problema:
L=1; K0=1; a=0;
%H1=3; H2=3; uE=0;
H1=0; H2=0; uE=0;
% Q =@(x) (10*exp(-100*(x-0.5).^2));
% Q=@(x) (100*(x<0.7));
 Q=@(x) 100*(1-x);

% parametros del nodos de resolucion
N=10; % cantidad de subintervalos

%% Armado de la matriz
h = L/N;
unos = ones(N+2,1);
diagonales = [-1*unos 2*unos -1*unos];
matriz = spdiags(diagonales, [-1 0 1], N+2, N+2);
matriz(1,[1:2])=[1 0];% Nodo izq condicion DIRICLET
matriz(N+2, [N:N+2]) = [-1 2*h*H1/K0 1];

% Armado del lado derecho

x = linspace(0,L,N+1)';
F =[a; h^2/K0*Q(x(2:N+1)); 2*h*H2*uE/K0];

% Resolucion

U = matriz\F;

U(N+2) = []; % otra forma: U = U(1:N+1);
% Otra forma
%
Uinf = U(1:N+1);
figure(1)
plot(x,U,'*-')
grid on;
grid minor;
