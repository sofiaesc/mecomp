clf; clear all; clc;

# ----------------------------
# ------ EJERCICIO 1: --------

# ----------ítem b------------
inter = [0 2];
cb = [1, 100, -1; 2, 0, -1];
L = 50;
x = linspace(inter(1),inter(2),L);

k = 1;
rho = 1;
cp = 0;
cR = 1;
G = 0*x;
et = -1;
model = struct("k",k,"rho",rho,"cp",cp,"cR",cR,"G",G);

figure(1)
T = (100*e.^(-x).*(e.^(2*x)+e.^4))./(1+e.^4);    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = difFinitas_orden2(x, model, cb, et);
plot(x,Ta,'b*')
hold on
grid on
Tb = difFinitas_orden1(x, model, cb, et);
plot(x,Tb,'m--')


figure(2)
error(T,Tb)
hold on
grid on
error(T,Ta)
