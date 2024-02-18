% rutina para demostrar que la derivada tercera es de segundo orden

pkg load symbolic

syms x s
if 1
f = tanh(x^3/s);
else
f = 10*x^5;
end

d3fdx3 = diff(diff(diff(f,'x'),'x'),'x');

xp = 0.1; sp=0.1;
% valor analitico de la derivada tercera para sigma y x dado
d3fdx3_ex = double(subs(subs(d3fdx3,'s',sp),'x',xp))

% coeficientes
aa = [-3;-1;1;3];
dx = 0.1;

bfh = [];

for k=1:10

    if k>1,
        dx = dx/2;
    end

% evaluaciones de la funcion en los diferentes puntos
ff(1) = subs(subs(f,'s',sp),'x',xp+dx);
ff(2) = subs(subs(f,'s',sp),'x',xp-dx);
ff(3) = subs(subs(f,'s',sp),'x',xp+2*dx);
ff(4) = subs(subs(f,'s',sp),'x',xp-2*dx);

d3fdx3_num = double(simplify(ff*aa))/dx^3;

disp(' d3fdx3(ex) & (num) & error ')
error=abs(d3fdx3_num-d3fdx3_ex);
disp([num2str(d3fdx3_ex) '  &  '   num2str(d3fdx3_num) '  &  ' num2str(error)])

bfh = [bfh;[dx,error]]

end

figure,loglog(bfh(:,1),bfh(:,2))
xlabel(' \Delta x ','FontSize',14,'FontWeight','b');
ylabel(' error ','FontSize',14,'FontWeight','b');
title(' Estimacion numerica derivada 3ra a 2do orden ')
grid on
