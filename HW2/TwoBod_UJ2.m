function dx = TwoBod_UJ2(t, x)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector (orbit is Keplerian, doesn't matter) 
%   x = [6x1] state vector 
% 
% Outputs 
%   dx = [6x1] derivative of state vector 
% ------------------------------------------------------------------------

global mu J2 RE 

dx = zeros(6, 1);   % force column vector 

% dx1 = x4 
% dx2 = x5 
% dx3 = x6 
% dx4 = (-u/r^3) * x1
% dx5 = (-u/r^3) * x2 
% dx6 = (-u/r^3) * x3 

dx(1:3) = x(4:6); 
% r_norm  = sqrt( x(1)^2 + x(2)^2 + x(3)^2 ); 
r_norm  = norm(x(1:3)); 
dx(4:6) = ( - mu / r_norm^3 ) * x(1:3); 

% J2 stuff 
dUx =      -(3*J2*RE^2*mu*x*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)); 
dUy =      -(3*J2*RE^2*mu*y*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)); 
dUz =  -(3*J2*RE^2*mu*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)); 

dx(4:6) = dx(4:6) + [dUx; dUy; dUz]; 

end 