function dX = EOM(et, X)
% ------------------------------------------------------------------------
% Purpose: Generate EOM for satellite orbiting earth due to geopotential,
% lunisolar, SRP, and drag perturbations 
% 
% Inputs 
%   t   = [1x1] time (ET epoch) vector 
%   X   = [7x1] state vector in ECI frame (inertial) 
% 
% Outputs 
%   dX  = [7x1] derivative of state vector 
% ------------------------------------------------------------------------

global wE muE 
global A m p0 r0_drag H  

% force column vector. Check if X is numeric or symbolic 
if isnumeric(X)
    dX = zeros(7, 1);   
else 
    dX = sym(zeros(7,1)); 
end

% Set velocity and CD 
dX(1:3) = X(4:6);
CD = X(7); 

% accel due to point mass (not needed when geopotential gravity is present)
r       = norm(X(1:3)); 
% dX(4:6) = ( - muE / r^3 ) * X(1:3); 

% accel due to gravity
if isnumeric(X); g = fn.a_spherical(et, X); else g = fn.g_J2J3J4(X); end
% g = fn.a_spherical(et, X); 
% g = fn.g_J2J3J4(X); 
dX(4:6) = dX(4:6) + g; 

% accel due to lunisolar perturbation 
[a_sun, a_moon] = fn.lunisolar(et, X); 
dX(4:6) = dX(4:6) + a_sun + a_moon; 

% accel due to SRP 
a_srp   = fn.a_SRP(et, X); 
dX(4:6) = dX(4:6) + a_srp; 

% CHECK UNITS. USE KM !!!
% accel due to drag 
pA      = p0 * exp( -(r - r0_drag)/H ); 
VA      = X(4:6) - cross( [0; 0; wE], X(1:3) ); 
VAnorm  = norm(VA); 
a_drag  = - 1/2 * CD * A/m * pA * VAnorm * VA;
% a_drag  = a_drag / 26.5;  % Correct to match Jah's Amat??? 
dX(4:6) = dX(4:6) + a_drag; 

end 