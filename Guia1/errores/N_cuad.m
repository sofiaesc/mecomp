# Con este script obtenemos los coeficientes para cada función de ajuste:
#
# K =

x1 = 2;
x1_2 = 2.5;
x2 = 3;

K = [x1.^2 x1 1;
     x1_2.^2 x1_2 1;
     x2.^2 x2 1];
coef = zeros(3,3);

for i=1:3
  F = zeros(3,1);
  F(i) = 1;
  coef_loc = K\F;
  coef(i,1:3) = coef_loc;
 endfor

 display(coef)

# Verificación:

# Suma de las a debe dar 0; suma de las b debe dar 0; suma de las c debe dar 1:
suma = sum(coef,1)  # Debo obtener 0, 0, 1

# N1 evaluado en x1 debe dar 1, en x1_2 o x2 debe dar 0:
disp('')
disp('N1')
coef_N1 = coef(1,:);
N1_x1 = polyval(coef_N1,x1)
N1_x1_2 = polyval(coef_N1,x1_2)
N1_x2 = polyval(coef_N1,x2)

# N1 evaluado en x1 debe dar 1, en x1_2 o x2 debe dar 0:
disp('')
disp('N2')
coef_N2 = coef(2,:);
N2_x1 = polyval(coef_N2,x1)
N2_x1_2 = polyval(coef_N2,x1_2)
N2_x2 = polyval(coef_N2,x2)

# N1 evaluado en x1 debe dar 1, en x1_2 o x2 debe dar 0:
disp('')
disp('N3')
coef_N3 = coef(3,:);
N3_x1 = polyval(coef_N3,x1)
N3_x1_2 = polyval(coef_N3,x1_2)
N3_x2 = polyval(coef_N3,x2)
