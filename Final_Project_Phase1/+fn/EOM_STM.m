function dX = EOM_STM(et, X, Amat_fn)
% ------------------------------------------------------------------------
% Purpose: Generate EOM for satellite orbiting earth due to gravity (EGM 96),
% lunisolar perturbations, SRP, and drag 
% 
% Inputs 
%   t   = [7x1] time (ET epoch) vector 
%   rv  = [49x1] state vector in ECI frame (inertial) 
% 
% Outputs 
%   drv = [49x1] derivative of state vector 
% ------------------------------------------------------------------------

% force column vector 
dX = zeros(7+49, 1);   

% EOM 
dX(1:7) = fn.EOM(et, X); 

% STM stuff 
Amat = Amat_fn(X(1), X(2), X(3), X(4), X(5), X(6), X(7)); 
STM  = X(8:7+49); 
STM  = reshape(STM, [7 7]); 
dSTM = Amat*STM; 
dSTM = reshape(dSTM, [49,1]); 

dX(8:7+49) = dSTM; 

end 



