function [T] = difFinitas_irregular(xnode, model, cb, ct)

  # ------- Constantes del modelo ------- #
  k = model.k;
  cR = model.cR;
  rho = model.rho;
  cp = model.cp;

  # ------- Armado del sistema de ecuaciones ------- #

  #  XNODE: Puntos de la malla
  #  Ejemplo: [0.0, 0.25, 0.3, 0.8, 0.9, 1]
  #  Resuelvo como si fuera irregular siempre, si h+ = h-
  #  se arma como regular porque vuelve al stencil clásico.

  K = zeros(length(xnode),length(xnode));
  for i = 2:(length(xnode)-1)
      hp = xnode(i+1) - xnode(i);   # h+
      hm = xnode(i) - xnode(i-1);   # h-

      a = 2/(hp*(hp+hm));
      b = -2/(hp*hm);
      c = 2/(hm*(hp+hm));

      A = (k*a);
      B = -(k*b+cr);
      C = (k*c);

      K(i,i-1:i+1) = [A,B,C];
      F(i) = model.G(i);
  endfor

  # ------- Condiciones de borde ------- #
  # Agrego primera y última fila según las condiciones de borde.
  # Usamos nodo ficticio, como está centrado lo hacemos como regular.

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

  # ------- Resolución del sistema ------- #
  T = K\F;

