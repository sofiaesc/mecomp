% Funcion para calcular los coeficientes asociados a un operador en diferencias
% dada la expresion dkT/dxk(i) = sum(aj*T(i+-j))

% puntos es un vector que indica que puntos alrededor de T(i) queremos usar
% y k es el orden de la derivada a aproximar, suponiendo una MALLA UNIFORME
% por ejemplo d3T/dx3 = aT(i+4)+bT(i+1)+cT(i-1)+dT(i-3) se calcula como:
% [coef] = taylorDF([4 1 -1 -3],3)
function [coef] = taylorDFsyms(puntos,k)
    if (length(puntos) <= k)
        sprintf('Faltan puntos para aproximar el operador de orden %i',k)
        coef = NaN;
    else
        syms h;
        N=length(puntos);
        M=zeros(N)*h;
        M(1,:)=1;
        for i=2:N
            M(i,:)=((puntos*h).^(i-1))./factorial(i-1);
        end
        b=zeros(N,1);
        b(k+1)=1;
        coef=M\b;
    end
end
