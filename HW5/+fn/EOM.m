function dX = EOM(et, X)
% ------------------------------------------------------------------------
% Purpose: Generate EOM for satellite orbiting earth due to gravity (EGM 96),
% lunisolar perturbations, SRP, and drag 
% 
% Inputs 
%   t   = [1x1] time (ET epoch) vector 
%   X   = [7x1] state vector in ECI frame (inertial) 
% 
% Outputs 
%   dX  = [7x1] derivative of state vector 
% ------------------------------------------------------------------------

global muE wE 
global A m p0 r0_drag H  

% force column vector 
dX = zeros(7, 1);   
dX(1:3) = X(4:6); 
CD = X(7); 

% accel due to point mass 
r       = norm(X(1:3)); 
dX(4:6) = ( - muE / r^3 ) * X(1:3); 

% accel due to gravity - DONE 
g_J2J3J4 = fn.g_J2J3J4(X); 
g_ECI    = fn.a_spherical(et, X); 
dX(4:6)  = dX(4:6) + g_ECI; 

% accel due to lunisolar perturbation - DONE 
[a_sun, a_moon, X_ESun, ~] = fn.lunisolar(et, X); 
dX(4:6) = dX(4:6) + a_sun + a_moon; 

% accel due to SRP 
a_srp = fn.a_SRP(X, X_ESun); 
dX(4:6) = dX(4:6) + a_srp; 

% accel due to drag 
pA      = p0 * exp( -(r - r0_drag)/H ); 
VA      = [X(4) + wE * X(2); ... 
           X(5) - wE * X(1); ... 
           X(6)]; 
VAnorm  = norm(VA); 
a_drag  = - 1/2 * CD * A/m * pA * VAnorm * VA;  
dX(4:6) = dX(4:6) + a_drag; 

end 



%% functions 

function rv = spice_state(epoch, target, frame, abcorr, observer) 

    rv = zeros(length(epoch), 6); 
    
    for i = 1:length(epoch) 

        %  Look-up the state for the defined parameters.
        starg   = mice_spkezr( target, epoch(i), frame, abcorr, observer);
        rv(i,:) = starg.state(1:6); 
        
    end 

end 


