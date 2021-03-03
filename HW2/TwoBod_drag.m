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

% dx1 = x4 
% dx2 = x5 
% dx3 = x6 
% dx4 = (-u/r^3) * x1
% dx5 = (-u/r^3) * x2 
% dx6 = (-u/r^3) * x3 

dx(1:3) = x(4:6); 
% r_norm  = sqrt( x(1)^2 + x(2)^2 + x(3)^2 ); 
rnorm  = norm(x(1:3)); 
dx(4:6) = ( - mu / rnorm^3 ) * x(1:3); 

% debug drag stuff 
rv = x; 
x = rv(1); 
y = rv(2); 
z = rv(3); 
xd = rv(4); 
yd = rv(5); 
zd = rv(6); 


% drag stuff 
pA = p0 * exp( -(rnorm - r0_drag)/H ); 
% V_A = [dx(4) + dtheta * dx(2); dx(5) - dtheta * dx(1); dx(6)]; 
V_A = [xd + dtheta * y; yd - dtheta * x; zd]; 
% VA = sqrt( ( dx(4) + dtheta * dx(2) )^2 + ( dx(5) - dtheta * dx(1) )^2 + dx(6)^2 ); 
VA = sqrt( ( xd + dtheta*y )^2 + ( yd - dtheta*x )^2 + zd^2 ); 

acc_drag = - 1/2 * CD * A/m * pA * VA * V_A; 
dx(4:6) = dx(4:6) + acc_drag; 

end 