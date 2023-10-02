# para saber orden del error hago log(RMSE) y busco una recta que aproxime
# mejor los puntos (polyval en octave). de ahi saco el p y la ordenada al
# origen de la recta. podría no poner nada en el eje x pero es el log de los dx

# script que vaya graficando el error con distintos puntos (8, 16, ...)
# con respecto a la analitica para ir viendo

function error(yd,yp)
  rmse = sqrt(mean((yd - yp).^2));
  p = polyfit(yd,yp);

  plot(p)
endfunction

