function [T] = difFinitas_orden1(xnode, model, cb, et)

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
  # Usamos primer orden de aproximación (sin nodo ficticio).
  b = zeros(length(xnode),1);
  b(1:end) = (G(1:end)./k)*(dx^2);

  # BORDE A:
  switch(cb(1,1))
    case 1 # Dirichlet
      Td=cb(1,2);
      K(1,1) = 1;
      b(1) = Td;
    case 2 # Neumann
      q=cb(1,2);
      K(1,1) = -1;
      K(1,2) = 1;
      b(1) = -(dx*q)./k;
    case 3 # Robin
      h = cb(1,2);
      uE = cb(1,3);
      K(1,1) = -k./dx + h;
      K(1,2) = k./dx;
      b(1) = h*uE;
  endswitch

  # BORDE B:
  switch(cb(2,1))
    case 1 # Dirichlet
      Td = cb(2,2);
      K(end,end) = 1;
      b(end) = Td;
    case 2 # Neumann
      q = cb(2,2);
      K(end,end) = 1;
      K(end,end-1) = -1;
      b(end) = -(dx*q)./k;
    case 3 # Robin
      h = cb(2,2);
      uE = cb(2,3);
      K(end,end) = k./dx + h;
      K(end,end-1) = -k./dx;
      b(end) = h*uE;
  endswitch

  # ------- Resolución del sistema ------- #
  T = K\b;

