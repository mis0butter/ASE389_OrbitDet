function dx = TwoBod_UJ2(t, x)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector  
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


% drag stuff 
pA = p0 * exp( -(rnorm - r0_drag)/H ); 
V_A = [dx(4) + dtheta * dx(2); 
    dx(5) - dtheta * dx(1); 
    dx(6)]; 
VA = sqrt( ( dx(4) + dtheta * dx(2) )^2 + ( dx(5) - dtheta * dx(1) )^2 + dx(6)^2 ); 

ddx_drag = - 1/2 * CD * A/m * pA * VA * V_A; 
dx(4:6) = dx(4:6) + ddx_drag; 

end 