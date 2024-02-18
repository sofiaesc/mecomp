function [T] = difFinitas(xnode, model, cb, et, Tant = zeros(N,1);)
  % inicializacion de variables %
  N = length(xnode);
  K = zeros(N,N);
  F = zeros(N,1);
  t = 0;

  % parametros del modelo %
  rho_cp = model.rho*model.cp;
  k = model.k;
  cr = model.c;

  % ecuaciones para nodos interiores %
  for i=2:N-1
    dx_izq = xnode(i)-xnode(i-1);
    dx_der = xnode(i+1)-xnode(i);

    a = 2/(dx_izq*(dx_der+dx_izq));
    b = 2/(dx_der*dx_izq);
    c = 2/(dx_der*(dx_der+dx_izq));

    A = (k*a);
    B = -(k*b+cr);
    C = (k*c);

    K(i,i-1:i+1) = [A B C];
    F(i) = model.G(i);
  endfor

  % Cond. de borde en x=L1 %
  dx = xnode(2)-xnode(1);
  d2x = dx*dx;
  if cb(1,1) == 2
    # Cond neumann
    q = cb(1,2);
    B = -((2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)-2*q/dx);
  elseif cb(1,1) == 3
    # Cond robin
    h = cb(1,2);
    Tinf = cb(1,3);
    B = -((2*h)/dx+(2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)+(2*h*Tinf)/dx);
  endif

  % Cond. de borde en x=L2 %
  dx = xnode(N)-xnode(N-1);
  d2x = dx*dx;
  if cb(2,1) == 2
    # Cond neumann
    q = cb(2,2);
    A = (2*k)/d2x;
    B = -((2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)-2*q/dx);
  elseif cb(2,1) == 3
    # Cond robin
    h = cb(2,2);
    Tinf = cb(2,3);
    A = (2*k)/d2x;
    B = -((2*h)/dx+(2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)+(2*h*Tinf)/dx);
  endif

  if et(1) == 0 # estado estacionario
    F *= -1;
    [K, F] = dirchlet(K, F, cb);
    T = inv(K)*F;
    plot(xnode, T, 'b-x')
    disp('Estado estacionario');

  elseif et(1) == 1 # esquema explicito
    T = Tant;
    maxIt = et(2); tol = et(3);
    dx = max(diff(xnode));
    alpha = k/rho_cp;
    dt = 0.9*((0.5*(dx*dx))/alpha); # 0.9*dt_critico

    [K, F] = dirchlet(K, F, cb);


    for i=1:maxIt
      Tact = (dt/rho_cp)*(K*Tant+F)+Tant;

      # en cada iteracion corregimos/fijamos la/las
      # condiciones dirichlet
      if cb(1,1) == 1
        Tact(1) = cb(1,2);
      endif
      if cb(2,1) == 1
        Tact(N) = cb(2,2);
      endif

      T = [T Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      plot(xnode, Tact, '*-b')
      title(sprintf('nit: %d - error: %e',i,err));
      pause(0.05)

      if (err<tol)
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;
    endfor
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario');
  else # esquema implicito
    T = Tant;
    maxIt = et(2); tol = et(3);
    dx = max(diff(xnode));
    dt = et(4);

    Kimp = (rho_cp/dt)*eye(N,N) - K;
    for i=1:maxIt
      Fimp = F + (rho_cp/dt)*Tant;
      Tact = inv(Kimp)*Fimp;

      # en cada iteracion corregimos/fijamos la/las
      # condiciones dirichlet
      if cb(1,1) == 1
        Tact(1) = cb(1,2);
      endif
      if cb(2,1) == 1
        Tact(N) = cb(2,2);
      endif

      T = [T Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      plot(xnode, Tact, '*-b')
      title(sprintf('nit: %d - error: %e',i,err));
      pause(0.05)

      if (err<tol)
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;
    endfor
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario');

  endif

endfunction
