%% (d^2T/dx^2)+100*((x-3)^2) = 0                                      %%
%% dT/dx(0) = -2; T(5) = 0                                            %%
%% Sol. exacta: T(x) = (-25*(x^4)+300*(x^3)-1350*(x^2)+1906*x+2345)/3 %%
clf; clear all; clc;
pkg load symbolic
warning('off', 'all');

syms x;
T_ex = (-25*(x^4)+300*(x^3)-1350*(x^2)+1906*x+2345)/3;
x_num = 1:0.1:5;

%% FORMUL. FUERTE %%
M = 6; % cantidad de funciones de prueba
N = x.^(1:M);
dN = diff(N,1);
d2N = diff(N,2);

K = int(N'*d2N,1,5)+subs(N'*dN,x,1)-subs(N'*N,x,5);
F = -100*int(N*((x-3)^2),1,5)-subs(2*N,x,5);
a = double(K)\double(F)';
T_f = N*a;

subplot(1,2,1);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_f,x,x_num)),'rx-');
title('Formulacion fuerte; M=6');

%% FORMUL. DEBIL %%
M = 6; % cantidad de funciones de prueba
N = x.^(1:M);
dN = diff(N,1);

K = int(dN'*dN,1,5)-subs(N'*dN,x,5)+subs(N'*N,x,5);
F = 100*int(N*((x-3)^2),1,5)+subs(2*N,x,1);
a = double(K)\double(F)';
T_d = N*a;

subplot(1,2,2);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_d,x,x_num)),'go-');
title('Formulacion debil; M=6');
