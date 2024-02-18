function [K,F] = fdm2d_dirichlet(K,F,DIR)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Dirichlet.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * DIR: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: número de nodo.
%   - Columna 2: valor en ese nodo (escalar)

% Salida:
% * K: matriz del sistema (difusión + reacción) luego de realizar las simplificaciones
%   que surgen de aplicar la condición de borde Dirichlet.
% * F: vector de flujo térmico luego de realizar las simplificaciones que surgen de
%   aplicar la condición de borde Dirichlet.
% ----------------------------------------------------------------------

%Cantidad de nodos de condicion Dirichlet
P = size(DIR);

for i=1:P
    %Obtenemos el nodo y su valor de condicion
    nodo = DIR(i,1);
    valor = DIR(i,2);

    %Ajustamos la matriz K
    K(nodo,:)=0;
    K(nodo,nodo)=1;

    %Ajustamos el vector F
    F(nodo)=valor;

end

end

