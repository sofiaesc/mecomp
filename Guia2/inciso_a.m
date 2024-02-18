%% 2*(d^2T/dx^2)+100 = 0                 %%
%% T(0) = 10; T(1) = 50                  %%
%% Sol. exacta: T(x) = -25*(x^2)+65*x+10 %%

pkg load symbolic
warning('off', 'all');

syms x;
T_ex = -25*(x^2)+65*x+10;
x_num = 0:0.05:1;

%% CON PSI %%
psi = 40*x+10;
M = 2; % cantidad de funciones de prueba
N = sin((1:M)*sym(pi)*x);
d2N = diff(N,2);

K = 2*int(N'*d2N,0,1);
F = -100*int(N,0,1);
a = double(K)\double(F)';
T_ap1 = psi+N*a;

subplot(1,2,1);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_ap1,x,x_num)),'rx-');
title('Psi y Nm = sin(m*pi*x); M=2');
ylim([0, 55]);

%% SIN PSI %%
M = 8; % cantidad de funciones de prueba
N = exp((1:M)*x);
d2N = diff(N,2);

K = 2*int(N'*d2N,0,1)-subs(N'*N,x,0)-subs(N'*N,x,1);
F = -100*int(N,0,1)-10*subs(N,x,0)-50*subs(N,x,1);
a = double(K)\double(F)';
T_ap2 = N*a;

subplot(1,2,2);
plot(x_num,double(subs(T_ex,x,x_num)));
hold on;
plot(x_num,double(subs(T_ap2,x,x_num)),'go-');
title('Nm = e^{m*x}; M=8');
ylim([0, 55]);
