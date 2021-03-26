function dx = TwoBod_4states(t, x)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector (orbit is Keplerian, doesn't matter) 
%   x = [4x1] state vector 
% 
% Outputs 
%   dx = [4x1] derivative of state vector 
% ------------------------------------------------------------------------

global mu 

dx = zeros(4, 1);   % force column vector 

% dx1 = x3 
% dx2 = x4 
% dx3 = (-u/r^3) * x1
% dx4 = (-u/r^3) * x2

dx(1:2) = x(3:4); 
r_norm  = norm(x(1:2)); 
dx(3:4) = ( - mu / r_norm^3 ) * x(1:2); 

end 