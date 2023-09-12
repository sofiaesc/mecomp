function [T] = difFinitas_orden2(xnode, model, cb, et)

  # ------- Constantes del modelo ------- #
  k = model.k;
  cR = model.cR;
  rho = model.rho;
  cp = model.cp;
  G = model.G;

  # ------- Armado del sistema de ecuaciones ------- #

  # FOR PARA ARMAR LA MATRIZ (malla regular):
  dx = xnode(2)-xnode(1);
  K = zeros(length(xnode),length(xnode));
  for i = 2:(length(xnode)-1)
    K(i,i-1:i+1) = [-1 2+((dx.^2*cR)./k) -1];
  endfor

  # ------- Condiciones de borde ------- #
  # Agrego primera y última fila según las condiciones de borde.
  # Usamos nodo ficticio (2do orden de aproximación).
  b = zeros(length(xnode),1);
  b(1:end) = (G(1:end)./k)*(dx^2);

  # BORDE A:
  switch(cb(1,1))
    case 1  # Dirichlet
      Td=cb(1,2);
      K(1,1) = 1;
      b(1) = Td;
    case 2  # Neumann
      q=cb(1,2);
      K(1,1) = 2 + (dx.^2*cR)./k;
      K(1,2) = -2;
      b(1) += -(2*dx*q)./k;
    case 3  # Robin
      h = cb(1,2);
      uE = cb(1,3);
      K(1,1) = 2 + (dx.^2*cR)./k + (2*dx*h)./k;
      K(1,2) = -2;
      b(1) += ((2*dx*h)./k).*uE;
  endswitch

  # BORDE B:
  switch(cb(2,1))
    case 1  # Dirichlet
      Td = cb(2,2);
      K(end,end) = 1;
      b(end) = Td;
    case 2  # Neumann
      q = cb(2,2);
      K(end,end) = 2 + (dx.^2*cR)./k;
      K(end,end-1) = -2;
      b(end) += -(2*dx*q)./k;
    case 3  # Robin
      h = cb(2,2);
      uE = cb(2,3);
      K(end,end) = 2 + (dx.^2*cR)./k + (2*dx*h)./k;
      K(end,end-1) = -2;
      b(end) += ((2*dx*h)./k).*uE;
  endswitch

  # ------- Condición temporal ------- #
  # Consideramos siempre la condición nula: T(x,0) = 0.

  #switch(ct)
   # case -1:  # Estacionario
   # case 1:   # Forward Euler
   # case 2:   # Backward Euler
   # case 3:   # Crank-Nicholson

  # for (
      # T = K\b
      # El K no cambia, el b sí cambia y por eso está en un for.

  # ------- Resolución del sistema ------- #
  T = K\b;

