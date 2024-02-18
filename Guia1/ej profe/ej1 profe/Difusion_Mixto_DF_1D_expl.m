% Difusion_Mixto_DF_1D_expl.m
%
% Resuelve la ecuacion de Difusion
% con condiciones de Dirichlet en x=0 y Mixtas en x=L.
% Metodo: Diferencias finitas en espacio y Euler EXPLICITO en tiempo.
%
% c rho u_t(x,t) - K0 u_xx(x,t) = Q(x,t), 0 < x < L, 0 < t < T,
% u(0,t) = a(t), K0 u_x(L,t) + H1 u(L,t) = H2 uE(t), 0 < t < T,
% u(x,0) = u0(x), 0 < x < L.
%
% Calculamos la solucion estacionaria
Poisson_DR

%%Parametros del problema
L = 1; % longitud del dominio
T = 0.3; % PERIODO DE CALCULO
K0 = 0.5; % conductividad
c = 1; % Calor especifico
rho = 1; % densidad
% FUENTES
%Q = @(x,t)(10*exp(-100*(x-.5).^2)*(t<0.25));
Q = @(x,t) (100*(1-x))*(t<10);

% CONDICION DE FRONTERA
a = @(t)(0);
H1 = 0; H2 = 0; uE = @(t)(0);
u0 = @(x)(zeros(size(x)));
%u0 = @(x) (5*x.*(2-x));


%%Constantes y funciones utiles
k = K0/(c*rho); % difusividad
q = @(x,t) Q(x,t) / (c*rho); % VALOR INTESIVO DE LA FUENTE

%%Parametros del metodo de resolucion
N = 10; % CANTIDAD DE SUBINTERVALOS A DIVIDIR EL DOMINIO
h = L/N; % PASO

% ESTABILIDAD
lambda = 0.25; % la mitad de lo necesario para tener estabilidad
%lambda = 1 % No es estable
deltat = lambda*h^2/k
%==============================
%%Condicion inicial
X = linspace(0,L,N+1)';
U = u0(X);
% le damos valor a la U en el nodo ficticio (t=0)
U(N+2) = (H2*uE(0)-H1*U(N+1))*2*h/K0 + U(N);

% CONDICION INICIAL
ejes = [0 L 0 1];
figure(1); plot(X,U(1:N+1),'*-');
title('t = 0')
grid on; grid minor
pause(0.01);

%%Resolucion
t = deltat;
while (t < T)
% calculamos U en los nodos "interiores"
U(2:N+1) = deltat*q(X(2:N+1), t) ...
+ lambda*U(1:N) + (1-2*lambda)*U(2:N+1) + lambda*U(3:N+2);
% le damos valor a U en el extremo izquierdo
U(1) = a(t);
% le damos valor a la U en el nodo ficticio
U(N+2) = (H2*uE(t)-H1*U(N+1))*2*h/K0 + U(N);
figure(1);
plot(X,U(1:N+1),'*-');
title(sprintf('t = %5.3f',t));
grid on; grid minor
%axis(ejes);
pause(0.01);

if max(norm(U(1:N+1)-Uinf))<=0.001
  break
  elseif
  t = t +deltat
 endif
end

max(norm(U(1:N+1)-Uinf))

figure(2)
  plot(X,U(1:N+1),'*-b',X,Uinf,'r-o');
  title(sprintf('t = %5.3f',t));
  grid on, grid minor;



