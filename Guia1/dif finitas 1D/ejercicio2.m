clf; clear all; clc;

# ----------------------------
# ------ EJERCICIO 1: --------

# ----------ítem a------------
inter = [0 1];
cb = [1, 10, -1; 1, 50, -1];
L = 50;
x = linspace(inter(1),inter(2),L);

k = 2;
rho = 1;
cp = 0;
cR = 0;
G = 100*x.^0;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(1)
T = -25*x.^2 + 65*x + 10;    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')

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
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(2)
T = (100*e.^(-x).*(e.^(2*x)+e.^4))./(1+e.^4);    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')


# ----------ítem c------------
inter = [1 5];
cb = [2,2,-1; 1,0,-1];
L = 50;
x = linspace(inter(1),inter(2),L);

rho = 1;
cp = 0;
k = 1;
cR = 0;
et = [0 -1 -1 -1];
G = 100*(x-3).^2;
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(3)
T = 1/3*(-25*x.^4+300*x.^3-1350*x.^2+1906*x+2345);    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')


# ----------ítem d------------
inter = [0 1];
cb = [1,10,-1; 3,0.2,50];
L = 50;
x = linspace(inter(1),inter(2),L);

rho = 1;
cp = 0;
k = 1;
cR = 1;
G = 50*x.^0;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(4)
T = -36.6897*(e.^-x)-3.3103*e.^x+50;    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')


# ----------ítem e------------
inter = [5 10];
cb = [3,2,100; 1,50,-1];
L = 50;
x = linspace(inter(1),inter(2),L);

rho = 1;
cp = 1;
k = 2;
cR = 0;
G = x.^3;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(5)
T = (-1/40)*x.^5 + (1225/3)*x - 4600/3;    # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')


# ----------ítem f------------
inter = [0 1];
cb = [1,0,-1; 3,2,10];
L = 50;
x = linspace(inter(1),inter(2),L);

rho = 1;
cp = 2;
k = 2;
cR = 2;
G = 75*x.^0;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(6)
T = -5/4*e.^(-(x+1)).*(e.^x-1).*(11*e.^x+11-30*e);   # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')

# ----------ítem g------------
inter = [0 1];
cb = [1,50,-1; 2,5,-1];
L = 50;
x = linspace(inter(1),inter(2),L);

rho = 1;
cp = 1;
k = 2;
cR = -2;
G = 0*x;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(7)
T = 73.2433*sin(x) + 50*cos(x);  # Solución analítica
plot(x,T,'k')
hold on
grid on

Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')
