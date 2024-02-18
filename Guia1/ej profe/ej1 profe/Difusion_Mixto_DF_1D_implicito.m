% Difusion_Mixto_DF_1D_implicito.m
%
% Resuelve la ecuacion de Difusion
% con condiciones de Dirichlet en x=0 y Robin en x=L.
% Metodo: Diferencias finitas en espacio y Euler IMPLICITO en tiempo.
%
% u_t(x,t) - k u_xx(x,t) = f(x,t), 0 < x < L, t > 0,
% u(0,t) = a(t), K0 u'(L,t) + H1 u(L,t) = H2 uE(t), t > 0,
% u(x,0) = u0(x), 0 < x < L.

% Calculamos la solucion estacionaria
Poisson_DR

%%Parametros del problema
L = 1; T = 4;
K0 = 1; c=1; rho=1;
%f = @(x,t)(10*exp(-100*(x-.5).^2)*(t<1));
Q = @(x,t) (100*(1-x))*(t<10);
a = @(t)(0);
%H1 = 3; H2 = 3; uE = @(t)(0);
H1 = 0; H2 = 0; uE = @(t)(0);
u0 = @(x)(zeros(size(x)));
%u0 = @(x) (5*x.*(2-x));
% Constantes y Funciones Utiles
k=K0/(c*rho);
q = @(x,t) Q(x,t)/(c*rho);

%%Parametros del metodo de resolucion
N = 20;
h = L/N;
lambda = 20; % la mitad de lo necesario para tener estabilidad
deltat = lambda*h^2/k;

%%Condicion inicial
X = linspace(0,L,N+1)';
U = u0(X);
% Graficamos condicion Inicial
%ejes = [0 L 0 1];
figure(1); plot(X,U(1:N+1));
title('t = 0');
%axis(ejes);
pause

%%Armado de la matriz
unos = ones(N+2,1);
columnas = [-lambda*unos (1+2*lambda)*unos -lambda*unos];

matriz = spdiags(columnas, [-1 0 1], N+2, N+2);
matriz(1,1:3) = [1 0 0];
matriz(N+2, N:N+2) = [-1 2*h*H1/K0 1];

%%Resolucion
t=deltat;
while (t<T)
  % Ensamblado lado Derecho
  F = [a(t)
  U(2:N+1)+deltat*q(X(2:N+1),t)
  2*h*H2/K0*uE(t) ];
  % Calculamos la nueva U
  U = matriz\F;

  % Graficamos
  figure(1)
  plot(X,U(1:N+1),'*-b');
  title(sprintf('t = %5.3f',t));
  %axis(ejes);
  grid on, grid minor;
  pause(0.01)

if max(norm(U(1:N+1)-Uinf))<=0.01
  break
  elseif
  t = t +deltat
 endif
end

Error = max(norm(U(1:N+1)-Uinf))

figure(2)
  plot(X,U(1:N+1),'*-b',X,Uinf,'r-o');
  title(sprintf('t = %5.3f',t));
  grid on, grid minor;
