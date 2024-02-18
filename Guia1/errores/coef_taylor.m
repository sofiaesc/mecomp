pkg load symbolic
clc; clear all; clf;

phi = @(x) exp(-x.^2);
dphi = @(x) -2*x.*exp(-x.^2);
real = dphi(0.5);

h_ini = 0.1;
h = [h_ini h_ini/2 h_ini/4 h_ini/8];


K = [ 1 1 1 1;
      -2.1 -1.1 1.21 0;
     4.41 1.21 1.4641 0;
     -9.261 -1.331 1.772 0];
errores = [];

for i=1:length(h)
  nodos = [0.5-2.1*h(i) 0.5-1.1*h(i) 0.5+1.21*h(i) 0.5];
  val = phi(nodos)
  f = [0; 0; 1/2*h(i).^2; 0];
  coef = K\f;
  aprox = coef(1)*val(1) + coef(2)*val(2) + coef(3)*val(3) + (coef(1)+coef(2)+coef(3)+coef(4))*val(4);
  error_act = real - aprox;
  errores = [errores error_act];
endfor

figure(1)
errores
plot(log(h),log(errores))
poli = polyfit(log(h),log(errores),1);
p = poli(1)
