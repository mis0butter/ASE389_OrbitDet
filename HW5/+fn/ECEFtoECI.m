function [r_ECI] = ECEFtoECI(JD, r_ECEF)
% ------------------------------------------------------------------------ 
% Purpose: Convert ECF (ECEF/ITRF) position to ECI (ICRF) position 
% 
% Inputs: 
%   eop_data = IAU1980 EOP data (finals.all)
%   JD       = Julian Date (UTC) 
%   r_ECEF   = position in ECEF (Earth-centered Earth-fixed) frame
% 
% Outputs: 
%   r_ECI    = position in ECI (Earth-centered inertial) frame 
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

global eop_data 

% P  = precession of ECF wrt ECI 
P = fn.precession(JD); 
  
% N  = nutation of ECF wrt ECI  
[N, em, dpsi] = fn.nutation(JD); 

% S' = rotation of ECF about angular velocity vector 
[xp, yp, dT] = fn.iers_data(eop_data, JD); 
GSMT = 4.894961212823058751375704430 + dT * ... 
    ( 6.300388098984893552276513720 + dT * ... 
    ( 5.075209994113591478053805523e-15 - ... 
    -9.253097568194335640067190688e-24 * dT) ); 
 
aG = GSMT + dpsi * cos(em); 
Sp = [cos(aG), sin(aG), 0; -sin(aG), cos(aG), 0; 0, 0, 1 ]; 

% W  = offset of Earth's angular velocity vector wrt ECEF Z axis 
W = [1 0 xp; 0 1 -yp; -xp yp 1]; 

% ECI position calculation 
ECI_C_ECEF = (W * Sp * N * P)'; 
r_ECI      = fn.orthodcm(ECI_C_ECEF) * r_ECEF; 

end 