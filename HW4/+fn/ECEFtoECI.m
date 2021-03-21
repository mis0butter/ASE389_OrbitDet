function [r_ECI] = ECEFtoECI(JD, r_ECEF)
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
% Notes: 
%   Transformation to ECI from ECEF: 
%   ECI_DCM_ECEF = W * S' * N * P; 
%   W  = offset of Earth's angular velocity vector wrt ECEF Z axis 
%   S' = rotation of ECF about angular velocity vector 
%   N  = nutation of ECF wrt ECI 
%   P  = precession of ECF wrt ECI 
% ------------------------------------------------------------------------ 


%% P  = precession of ECF wrt ECI 

% time = number of centuries since J2000 as terrestrial time (TT) 
t = (JD - 2451545.0)./36525;

% precession angles??????
sigma = 2306.2181 * t + 0.30188 * t^2 + 0.017998 * t^3; 
theta = 2004.3109 * t - 0.42655 * t^2 - 0.041833 * t^3; 
z     = 2306.2181 * t + 1.09468 * t^2 + 0.018203 * t^3; 

% P row 1 coeffs 
p11 = cos(sigma)*cos(theta)*cos(z) - sin(sigma)*sin(z); 
p12 = -sin(sigma)*cos(theta)*cos(z) - cos(sigma)*sin(z); 
p13 = -sin(theta)*cos(z); 

% P row 2 coeffs 
p21 = cos(sigma)*cos(theta)*sin(z) + sin(sigma)*cos(z); 
p22 = -sin(sigma)*cos(theta)*sin(z) + cos(sigma)*cos(z); 
p23 = -sin(theta)*sin(z); 

% P row 3 coeffs 
p31 = cos(sigma)*cos(theta); 
p32 = -sin(sigma)*sin(theta); 
p33 = cos(theta); 

% P  = precession of ECF wrt ECI 
P = [ p11, p12, p13 ; 
      p21, p22, p23 ; 
      p31, p32, p33 ]; 
  
  
%% N  = nutation of ECF wrt ECI 




%% ECI position calculation 

ECI_C_ECEF = W * S' * N * P; 
r_ECI      = ECI_C_ECEF * r_ECEF; 

end 