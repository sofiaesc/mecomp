function [T] = tp_dif_finitas(xnode, model, cb, et, Tant = zeros(length(xnode),1))

  #--- Inicialización de variables ---#
  N = length(xnode);
  K = zeros(N,N);
  F = zeros(N,1);

  #--- Parámetros del modelo ---#
  rho_cp = model.rho*model.cp;
  k = model.k;
  cr = model.c;

  #--- Ecuaciones para nodos interiores ---#
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

  #--------------------------------------------#
  #--- Cond. de borde en x=L1 (extremo izq) ---#
  #--------------------------------------------#
  dx = xnode(2)-xnode(1);
  d2x = dx*dx;
  if cb(1,1) == 2 #--- Cond Neumann ---#
    q = cb(1,2);
    B = -((2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)-2*q/dx);
  elseif cb(1,1) == 3 #--- Cond Robin ---#
    h = cb(1,2);
    Tinf = cb(1,3);
    B = -((2*h)/dx+(2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)+(2*h*Tinf)/dx);
  endif

  #--------------------------------------------#
  #--- Cond. de borde en x=L2 (extremo der) ---#
  #--------------------------------------------#
  dx = xnode(N)-xnode(N-1);
  d2x = dx*dx;
  if cb(2,1) == 2 #--- Cond Neumann ---#
    q = cb(2,2);
    A = (2*k)/d2x;
    B = -((2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)-2*q/dx);
  elseif cb(2,1) == 3 #--- Cond Robin ---#
    h = cb(2,2);
    Tinf = cb(2,3);
    A = (2*k)/d2x;
    B = -((2*h)/dx+(2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)+(2*h*Tinf)/dx);
  endif

  if et(1) == 0
    #-----------------------------------------------#
    #--- Tratamiento para el estado estacionario ---#
    #-----------------------------------------------#
    F *= -1;
    [K, F] = dirchlet(K, F, cb);
    T = inv(K)*F;
    plot(xnode, T, 'b-x')
    disp('Estado estacionario');


  elseif et(1) == 1
    #---------------------------------------------#
    #--- Tratamiento para el esquema explicito ---#
    #---------------------------------------------#
    T = Tant;
    maxIt = et(2); tol = et(3);
    dx = max(diff(xnode));
    alpha = k/rho_cp;
    dt = et(4);
    [K, F] = dirchlet(K, F, cb);

    for i=1:maxIt

      Tact = (dt/rho_cp)*(K*Tant+F)+Tant;

      if cb(1,1) == 1       # Se pisa la cond Dirichlet en L1 si existe
        Tact(1) = cb(1,2);
      endif
      if cb(2,1) == 1       # Se pisa la cond Dirichlet en L2 si existe
        Tact(N) = cb(2,2);
      endif

      T = [T Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      if (err<tol)
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;
    endfor
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario');
  elseif et(1) == 2
    #---------------------------------------------#
    #--- Tratamiento para el esquema implicito ---#
    #---------------------------------------------#
    T = Tant;
    maxIt = et(2); tol = et(3);
    dx = max(diff(xnode));
    dt = et(4);
    Kimp = (rho_cp/dt)*eye(N,N) - K;

    for i=1:maxIt

      Fimp = F + (rho_cp/dt)*Tant;
      Tact = inv(Kimp)*Fimp;

      if cb(1,1) == 1       # Se pisa la cond Dirichlet en L1 si existe
        Tact(1) = cb(1,2);
      endif
      if cb(2,1) == 1       # Se pisa la cond Dirichlet en L2 si existe
        Tact(N) = cb(2,2);
      endif

      T = [T Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      if (err<tol)
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;

    endfor
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario');
  else
    #---------------------------------------------#
    #--- Tratamiento para el esquema semi (CN) ---#
    #---------------------------------------------#
    T = Tant;
    maxIt = et(2); tol = et(3);
    dx = max(diff(xnode));
    dt = et(4);
    Ksemi = (rho_cp/dt)*eye(N,N) - 1/2*K;

    for i=1:maxIt

      Fsemi = (1/2*K + (rho_cp/dt)*eye(N,N))*Tant + F;
      Tact = inv(Ksemi)*Fsemi;

      if cb(1,1) == 1       # Se pisa la cond Dirichlet en L1 si existe
        Tact(1) = cb(1,2);
      endif
      if cb(2,1) == 1       # Se pisa la cond Dirichlet en L2 si existe
        Tact(N) = cb(2,2);
      endif

      T = [T Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      if (err<tol)
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;

    endfor
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario')

  endif
endfunction
