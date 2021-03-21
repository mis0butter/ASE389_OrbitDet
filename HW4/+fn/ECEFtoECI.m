function [r_ECI] = ECEFtoECI(JD_TT, r_ECEF)
% ------------------------------------------------------------------------ 
% Purpose: Convert ECF (ECEF/ITRF) position to ECI (ICRF) position 
% 
% Inputs: 
%   JD      = Julan Date (UTC) 
%   r_ECEF  = position in ECEF (Earth-centered Earth-fixed) frame
% 
% Outputs: 
%   r_ECI   = position in ECI (Earth-centered inertial) frame 
% 
% References: 
%   Statistical Orbit Determination by Bob E. Schutz, George Born, and Tapley
% 
% Notes: 
%   Transformation to ECI from ECEF: 
%   ECI_DCM_ECEF = W * S' * N * P; 
%   W  = offset of Earth's angular velocity vector wrt ECEF Z axis 
%   S' = rotation of ECF about angular velocity vector 
%   N  = nutation of ECF wrt ECI 
%   P  = precession of ECF wrt ECI 
% ------------------------------------------------------------------------ 

% time = number of centuries since J2000 as terrestrial time (TT) 
t   = (JD_TT - 2451545.0)./36525;

%% P  = precession of ECF wrt ECI 

P = fn.precession(JD_TT); 
  
%% N  = nutation of ECF wrt ECI 

N = fn.nutation(JD_TT); 

%% S' = rotation of ECF about angular velocity vector 


%% W  = offset of Earth's angular velocity vector wrt ECEF Z axis 

W = [1 0 xp; 0 1 -yp; -xp yp 1]; 

%% ECI position calculation 

ECI_C_ECEF = W * S' * N * P; 
r_ECI      = ECI_C_ECEF * r_ECEF; 

end 