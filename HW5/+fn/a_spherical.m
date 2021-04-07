function g_ECI = a_spherical(et, X) 

global wE 

% accel due to 20x20 spherical harmonics gravity 

eop_data = load('finals_iau1980.txt'); 
JD_UTC   = cspice_et2utc(et, 'J', 10); 
JD_UTC   = str2num(extractAfter(JD_UTC, 'JD ')); 

% Convert r in ECI frame to ECEF frame. KM --> M
[r_ECEF] = fn.ECItoECEF(eop_data, JD_UTC, X(1:3));
r_ECEF   = r_ECEF*1000; 

% [g_ECEF] = fn.AccelHarmonic_ECEF(r_ECEF, 20, 20); 
[gx_ECEF, gy_ECEF, gz_ECEF] = gravitysphericalharmonic( [r_ECEF(1) r_ECEF(2) r_ECEF(3)], 'EGM96'); 
g_ECEF = [gx_ECEF, gy_ECEF, gz_ECEF]'; 

r_ECI  = X(1:3); 
v_ECI  = X(4:6); 
v_ECEF = v_ECI + cross([0, 0, wE], r_ECI)'; 

% Convert g accel in ECEF to ECI frame, accounting for rotating frames. M --> KM 
g_ECI = g_ECEF + [ -2*wE*v_ECEF(2); 2*wE*v_ECEF(1); 0 ] + [ -wE^2*r_ECI(1); -wE^2*r_ECI(2); 0 ]; 
% g_ECI = g_ECEF; 
g_ECI = g_ECI/1000; 

end 