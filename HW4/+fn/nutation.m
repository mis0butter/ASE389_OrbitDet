function N = nutation(JD_TT)

% time = number of centuries since J2000 as terrestrial time (TT) 
t      = (JD_TT - 2451545.0)./36525;
MJD_TT = JD_TT - 2400000.5; 
  
%% N  = nutation of ECF wrt ECI 

% em    = mean obliquity of the ecliptic 
% et    = true obliquity of the ecliptic 
% dpsi  = nutation in longitude 
% de    = nutation in obliquity 

% mean obliquity 
em = 84381.448 - 46.8150 * t - 0.00059 * t^2 + 0.001813 * t^3; 
em = em/3600 * pi/180; 

% nutation in longitude and obliquity ????? 
[dpsi, de] = nut_angles(JD_TT);

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

function [dpsi, deps] = nut_angles(JD_TT)

% JD_TT = Mjd_TT + 2400000.5; 

MJD_TT  = JD_TT - 2400000.5; 
T       = (MJD_TT - 51544.5)/36525;
rev     = 360*3600;  % arcsec/revolution

C = load('nut80.dat'); 

% Mean arguments of luni-solar motion

%   l   mean anomaly of the Moon
%   lp  mean anomaly of the Sun
%   F   mean argument of latitude
%   D   mean longitude elongation of the Moon from the Sun 
%   Om  mean longitude of the ascending node  
l  = mod (  485866.733 + (1325.0*rev +  715922.633)*T + 31.310*T^2 + 0.064*T^3, rev );
lp = mod ( 1287099.804 + (  99.0*rev + 1292581.224)*T -  0.577*T^2 - 0.012*T^3, rev );
F  = mod (  335778.877 + (1342.0*rev +  295263.137)*T - 13.257*T^2 + 0.011*T^3, rev );
D  = mod ( 1072261.307 + (1236.0*rev + 1105601.328)*T -  6.891*T^2 + 0.019*T^3, rev );
Om = mod (  450160.280 - (   5.0*rev +  482890.539)*T +  7.455*T^2 + 0.008*T^3, rev );
                           
% errata???

% Nutation in longitude and obliquity [rad]  
dpsi = 0;
deps = 0;

for i = 1:length(C)
  arg  =  ( C(i,1) * l + C(i,2) * lp + C(i,3) * F + C(i,4) * D + C(i,5) * Om ) / 3600 * pi/180;
  dpsi = dpsi + ( C(i,6) + C(i,7) * T ) * sin(arg);
  deps = deps + ( C(i,8) + C(i,9) * T ) * cos(arg);
end
    
dpsi = 1e-5 * dpsi / 3600 * pi/180;
deps = 1e-5 * deps / 3600 * pi/180;

end