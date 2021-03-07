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
global CD A m p0 r0_drag H dtheta 

dx = zeros(6, 1);   % force column vector 

dx(1:3) = x(4:6); 
% r_norm  = sqrt( x(1)^2 + x(2)^2 + x(3)^2 ); 
rnorm  = norm(x(1:3)); 
dx(4:6) = ( - mu / rnorm^3 ) * x(1:3); 

% J2 stuff 
dUx =      -(3*J2*RE^2*mu*x(1)*(x(1)^2 + x(2)^2 - 4*x(3)^2)) / (2*(x(1)^2 + x(2)^2 + x(3)^2)^(7/2)); 
dUy =      -(3*J2*RE^2*mu*x(2)*(x(1)^2 + x(2)^2 - 4*x(3)^2)) / (2*(x(1)^2 + x(2)^2 + x(3)^2)^(7/2)); 
dUz =  -(3*J2*RE^2*mu*x(3)*(3*x(1)^2 + 3*x(2)^2 - 2*x(3)^2)) / (2*(x(1)^2 + x(2)^2 + x(3)^2)^(7/2)); 

dx(4:6) = dx(4:6) + [dUx; dUy; dUz]; 

% drag stuff 
pA = p0 * exp( -(rnorm - r0_drag)/H ); 
V_A = [x(4) + dtheta * x(2); x(5) - dtheta * x(1); x(6)]; 
VA = norm(V_A); 
% VA = sqrt( ( dx(4) + dtheta * dx(2) )^2 + ( dx(5) - dtheta * dx(1) )^2 + dx(6)^2 ); 

ddx_drag = - 1/2 * CD * A/m * pA * VA * V_A; 
dx(4:6) = dx(4:6) + ddx_drag; 

end 