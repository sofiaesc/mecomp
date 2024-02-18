# Problema:
# 2*(d2T/dx2) - 2*T = 0

# Condiciones de borde:
# Dirichlett T = 50 para x = 0
# Nuemann q = 5 para x = 1

clf; clc; clear all;
pkg load symbolic
warning('off', 'all');

syms x;
T_ex = 73.2433*sin(x) + 50*cos(x);
x_num = 0:0.05:1;

%% FORMULACIÓN FUERTE CON PSI %%
psi = -5*x + 50;
M = 10; % cantidad de funciones de prueba
N = exp(x.^(1:M));
dN = diff(N,1);
d2N = diff(N,2);

K = int(dN'*dN,0,1) - int(N'*dN,0,1);
F = 5*subs(N,x,1);
a = double(K)\double(F)';
T_f = N*a;

figure(1)
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_f,x,x_num)),'rx-');
title('Formulacion fuerte; M=6');

