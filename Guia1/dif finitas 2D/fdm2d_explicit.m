function [PHI_vec, Q_vec] = fdm2d_explicit(K,F,xnode,neighb,model,dt)
    A = dt/(model.rho*model.cp);

    I = eye(model.nnodes,model.nnodes);
    
    PHI = model.PHI_n;
    PHI_n = model.PHI_n;
    PHI_vec = PHI;
    Q_vec = zeros(model.nnodes,2);
    
    for n = 1 : model.maxit
        PHI = A*F + (I - A*K)*PHI_n;

        err = norm(PHI-PHI_n,2)/norm(PHI,2);
        
        PHI_n = PHI;
        PHI_vec = [PHI_vec PHI];
        [Q] = fdm2d_flux(PHI,neighb,xnode,model.k);
        Q_vec = [Q_vec, Q];
        
        if err < model.tol
            disp('MÃ©todo terminado por tolerancia de error.');
            return;
        end
    end

    disp('MÃ©todo terminado por lÃ­mite de iteraciones.');
end