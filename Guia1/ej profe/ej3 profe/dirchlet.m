function [K, F] = dirchlet(K, F, cb)
  if cb(1,1) == 1
    K(1,1) = 1;
    F(1) = cb(1,2);
  endif 
    
  if cb(2,1) == 1
    N = length(F);
    K(N,N) = 1;
    F(N) = cb(2,2);
  endif  
endfunction