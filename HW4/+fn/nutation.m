function [N, em, dpsi] = nutation(JD)

% time = number of centuries since J2000 as terrestrial time (TT) 
t   = (JD - 2451545.0)./36525;
MJD = JD - 2400000.5; 
  
%% N  = nutation of ECF wrt ECI 

% em    = mean obliquity of the ecliptic 
% et    = true obliquity of the ecliptic 
% dpsi  = nutation in longitude 
% de    = nutation in obliquity 

% mean obliquity 
em = 84381.448 - 46.8150 * t - 0.00059 * t^2 + 0.001813 * t^3; 
em = em/3600 * pi/180; 

% nutation in longitude and obliquity ????? 
[dpsi, de] = fn.nut_angles(JD);

% true obliquity 
et = em + de; 

n11 = cos(dpsi); 
n12 = -cos(em) * sin(dpsi); 
n13 = -sin(em) * sin(dpsi); 

n21 = cos(et) * sin(dpsi);
n22 = cos(em) * cos(et) * cos(dpsi) + sin(em) * sin(et); 
n23 = sin(em) * cos(et) * cos(dpsi) - cos(em) * sin(et); 

n31 = sin(et) * sin(dpsi); 
n32 = cos(em) * sin(et) * cos(dpsi) - sin(em) * cos(et); 
n33 = sin(em) * sin(et) * cos(dpsi) + cos(em) * cos(et); 

N = [n11 n12 n13; n21 n22 n23; n31 n32 n33]; 

end

