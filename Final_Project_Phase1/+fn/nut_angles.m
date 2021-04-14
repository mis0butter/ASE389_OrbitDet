function [dpsi, deps] = nut_angles(JD)

% JD_TT = Mjd_TT + 2400000.5; 
MJD = JD - 2400000.5; 
T   = (MJD - 51544.5)/36525;
rev = 360;  % deg/revoluton 

C = load('nut80.dat'); 
C(:, 6:end) = C(:, 6:end)*10; 

% From errata: Delaunay parameters, mean arguments of luni-solar motion in deg
%   Mm = mean anomaly of the Moon
%   Ms = mean anomaly of the Sun
%   uM = mean argument of latitude
%   D  = mean longitude elongation of the Moon from the Sun 
%   Om = mean longitude of the ascending node  
Mm = 134.96298139 + ( 1325*rev + 198.8673981 )*T + 0.0086972*T^2 + 1.78e-5*T^3; 
Ms = 357.52772333 + ( 99*rev   + 359.0503400 )*T - 0.0001603*T^2 - 3.3e-6*T^3; 
uM = 93.27191028  + ( 1342*rev + 82.0175381  )*T - 0.0036825*T^2 + 3.1e-6*T^3; 
D  = 297.85036306 + ( 1236*rev + 307.1114800 )*T - 0.0019142*T^2 + 5.3e-6*T^3; 
Om = 125.04452222 - ( 5*rev    + 134.1362608 )*T + 0.0020708*T^2 + 2.2e-6*T^3; 

% Nutation in longitude and obliquity 
dpsi = 0;
deps = 0;

for i = 1:length(C)
  api  =  ( C(i,1) * Mm + C(i,2) * Ms + C(i,3) * uM + C(i,4) * D + C(i,5) * Om ) * pi/180;
  dpsi = dpsi + ( C(i,6) + C(i,7) * T ) * sin(api);
  deps = deps + ( C(i,8) + C(i,9) * T ) * cos(api);
end
    
dpsi = 1e-5 * dpsi / 3600 * pi/180;
deps = 1e-5 * deps / 3600 * pi/180;

end