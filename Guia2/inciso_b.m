%% 2*(d^2T/dx^2)+100 = 0                 %%
%% T(0) = 10; T(1) = 50                  %%
%% Sol. exacta: T(x) = -25*(x^2)+65*x+10 %%

clf; clear all; clc;
pkg load symbolic
warning('off', 'all');

syms x;
T_ex = (100*exp(-x).*(exp(2*x)+exp(4)))/(1+exp(4));
x_num = 0:0.1:2;

%% SIN PSI, FUERTE %%
M = 20; % cantidad de funciones de prueba
N = exp((1:M)*x);
dN = diff(N,1);
d2N = diff(N,2);

K = int(N'*d2N,0,2) - int(N'*N,0,2) + subs(N'*N,x,0) - subs(N'*dN,x,2);
F = 100*subs(N,x,0);
a = double(K)\double(F)';
T_ap2 = N*a;

subplot(1,3,1);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_ap2,x,x_num)),'go-');
title('Nm = e^{m*x}; M=8');

%% SIN PSI, DEBIL %%
M = 20; % cantidad de funciones de prueba
N = exp((1:M)*x);
dN = diff(N,1);
d2N = diff(N,2);

K = -int(dN'*dN,0,2) - int(N'*N,0,2) + subs(N'*dN,x,0) + subs(N'*N,x,0);
F = -100*subs(N,x,0);
a = double(K)\double(F)';
T_ap2 = N*a;

subplot(1,3,2);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_ap2,x,x_num)),'go-');
title('Nm = e^{m*x}; M=8');

%% CON PSI, FUERTE %%
M = 20; % cantidad de funciones de prueba
N = x.*(x-2).^(1:M)
psi = 100;
dN = diff(N,1);
d2N = diff(N,2);

K = int(dN'*d2N,0,2) - int(N'*N,0,2) + subs(N'*dN,x,0) + subs(N'*N,x,0);
F = int(N*psi,0,2);
a = double(K)\double(F)';
T_ap2 = psi + N*a;

subplot(1,3,3);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_ap2,x,x_num)),'go-');
title('Nm = e^{m*x}; M=8');

