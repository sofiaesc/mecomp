clf; clear all; clc;
pkg load symbolic
warning('off', 'all');

syms x;
T_ex = (-1/12)*x.^4+(sym(pi).^3-6)*x/3+100;
x_num = 1:0.1:3.14;

%% SIN PSI - FORMUL. FUERTE %%
M = 5; % cantidad de funciones de prueba
N = (x-100).^(1:M);
dN = diff(N,1);
d2N = diff(N,2);

K = int(N'*d2N,1,sym(pi))+subs(N'*N,x,0)+subs(N'*dN,x,sym(pi));
F = -int(N*(x.^2),1,sym(pi))+100*subs(N,x,0)-subs(N,x,sym(pi));
a = double(K)\double(F)'
T_f = N*a

figure(1)
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_f,x,x_num)),'rx-');
title('Formulacion fuerte; M=6');
