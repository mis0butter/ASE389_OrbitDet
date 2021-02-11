function [rv] = oe2rv(oe)
% ------------------------------------------------------------------------ 
% Purpose: Convert orbital elements and time past epoch to the classic 
% Cartesian position and velocity
% 
% Inputs: 
%   oe      = [6x1] or [1x6] orbital elements 
%   delta_t = t - t0 time interval 
%   mu      = Gravity * Mass (of Earth) constant 
% 
% Outputs: 
%   rv      = position and velocity state vector 
% ------------------------------------------------------------------------ 

% global delta_t 
global mu 

a       = oe(1); 
e       = oe(2); 
i       = oe(3); 
omega   = oe(4); 
LAN     = oe(5); 
% M0      = oe(6); 
nu      = oe(6); 

% nu is TRUE ANOMALY --> use Kepler's to calculate MEAN ANOMALY 
% E = 2*atan( sqrt( (1-e)/(1+e) ) * tan(nu/2) ); 
% M = M0 + sqrt( mu/a^3 ) * (delta_t); 
% E = keplerEq(M, e, eps); 
% E = kepler(M, e); 
% nu = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E/2) ); 

p = a * ( 1 - e^2 );            % intermediate variable 
r = p / ( 1 + e*cos(nu) );      % r_magnitude, polar coordinates 

% Perifocal position and velocity 
% r_pf    = zeros(3,1);  
% v_pf    = zeros(3,1); 
% r_pf(3) = 0; 
% v_pf(3) = 0; 
% r_pf(1) = r * cos(nu); 
% r_pf(2) = r * sin(nu); 
% v_pf(1) = -sqrt(mu/p) * sin(nu); 
% v_pf(2) = sqrt(mu/p) * (e + cos(nu)); 

r_pf = [ r * cos(nu); r * sin(nu); 0 ]; 
v_pf = [ -sqrt(mu/p) * sin(nu); sqrt(mu/p) * (e + cos(nu)); 0 ]; 

% Perifocal to ECI transformation, 3-1-3 rotation 
R11 = cos(LAN)*cos(omega) - sin(LAN)*sin(omega)*cos(i); 
R12 = -cos(LAN)*sin(omega) - sin(LAN)*cos(omega)*cos(i); 
R13 = sin(LAN)*sin(i); 

R21 = sin(LAN)*cos(omega) + cos(LAN)*sin(omega)*cos(i); 
R22 = -sin(LAN)*sin(omega) + cos(LAN)*cos(omega)*cos(i); 
R23 = -cos(LAN)*sin(i); 

R31 = sin(omega)*sin(i); 
R32 = cos(omega)*sin(i); 
R33 = cos(i); 

R = [R11 R12 R13; R21 R22 R23; R31 R32 R33]; 

% Transform perifocal to ECI frame 
r_vec = R * r_pf; 
v_vec = R * v_pf; 

% Position and state vector 
rv = [r_vec; v_vec]; 

end 

%% Kepler equation solvers 

function E = keplerEq(M,e,eps)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [rad]. 
   	En  = M;
	Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
	while ( abs(Ens-En) > eps )
		En = Ens;
		Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
    end
	E = Ens;
end

function E = kepler(M, e)
    f = @(E) E - e * sin(E) - M;
    E = fzero(f, M);  % <-- I would use M as the initial guess instead of 0
end