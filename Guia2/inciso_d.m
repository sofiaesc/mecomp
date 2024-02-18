clf; clear all; clc;
pkg load symbolic
warning('off', 'all');

syms x;
T_ex = -36.6897*exp(-x) - 3.3103*exp(x) + 50;
x_num = 0:0.05:1;

##%% FORMUL. FUERTE %%
M = 20; % cantidad de funciones de prueba
N = exp((1:M)*x);
dN = diff(N,1);
d2N = diff(N,2);

K = int(N'*d2N,0,1)-int(N'*N,0,1)+subs(N'*N,x,0)+subs(N'*dN,x,1)+0.2*subs(N'*N,x,1);
F = 50*int(N,0,1)+subs(10*N,x,0)+subs(10*N,x,1);
a = double(K)\double(F)';
T_f = N*a;

subplot(1,3,1);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_f,x,x_num)),'rx-');
title('Formulacion fuerte; M=6');

%% FORMUL. DEBIL %%
M = 8; % cantidad de funciones de prueba
N = exp((1:M)*x);
dN = diff(N,1);

K = int(dN'*dN,0,1)-int(N'*N,0,1)+subs(N'*N,x,0)+subs(N'*dN,x,0)-0.2*subs(N'*N,x,1);
F = 50*int(N,0,1)+subs(10*N,x,0)-subs(10*N,x,1);
a = double(K)\double(F)';
T_d = N*a

subplot(1,3,2);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_d,x,x_num)),'go-');
title('Formulacion debil; M=6');

%% CON PSI. FUERTE %%
M = 6; % cantidad de funciones de prueba
N = (x-10).^(1:M);
psi = 10;
dN = diff(N,1);
d2N = diff(N,2);

K = int(N'*d2N,0,1)-int(N'*N,0,1)-subs(N'*dN,x,1)-0.2*subs(N'*N,x,1);
F = 50*int(N,0,1)+int(10*N,0,1)-subs(10*N,x,1);
a = double(K)\double(F)';
T_f = N*a;

subplot(1,3,3);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_f,x,x_num)),'rx-');
title('Formulacion fuerte; M=6');
