
# ----------ítem a------------
inter = [-4 4];
cb = [1,20,-1; 3,2,10];
L = 21;
x = linspace(inter(1),inter(2),L);

k = 3;
rho = 0;
cp = 0;
cR = 0;
G = 100*x;
et = [0 -1 -1 -1];
model = struct("k",k,"rho",rho,"cp",cp,"c",cR,"G",G);

figure(2)
Ta = dif_finitas(x, model, cb, et);
plot(x,Ta,'m*')

x(11)

Ta(end)

