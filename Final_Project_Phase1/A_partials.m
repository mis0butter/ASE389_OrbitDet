%%  DERIVE A MATRIX 

global muE wE
global A m p0 r0_drag H  
X  = sym('X', [7 1]); 
dX = sym('dX', [7 1]); 
CD  = X(7); 

% force column vector 
dX(1:3) = X(4:6); 
et = et_t0; 

% accel due to point mass 
r       = sqrt( X(1)^2 + X(2)^2 + X(3)^2 ); 
dX(4:6) = ( - muE / r^3 ) * X(1:3); 

% accel due to gravity: J2, J3, J4 
g_J2J3J4 = fn.g_J2J3J4(X); 
dX(4:6) = dX(4:6) + g_J2J3J4; 

% accel due to lunisolar perturbation - DONE 
[a_sun, a_moon, X_ESun, X_EMoon] = fn.lunisolar(et, X); 
dX(4:6) = dX(4:6) + a_sun + a_moon; 

% accel due to SRP 
a_srp = fn.a_SRP(X, X_ESun); 
dX(4:6) = dX(4:6) + a_srp; 

% accel due to drag
pA       = p0 * exp( -(r - r0_drag)/H ); 
Vrel     = [X(4) + wE * X(2); ... 
            X(5) - wE * X(1); ... 
            X(6)]; 
Vnorm    = norm(Vrel); 
a_drag   = - 1/2 * CD * A/m * pA * Vnorm * Vrel;  
dX(4:6) = dX(4:6) + a_drag; 

% compute partials 
Amat = jacobian( dX, X );       
Amat; 

