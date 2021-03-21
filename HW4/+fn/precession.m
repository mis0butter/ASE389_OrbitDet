function P = precession(JD_TT) 

t   = (JD_TT - 2451545.0)./36525;

% precession angles ... in arcseconds 
zeta  = 2306.2181 * t + 0.30188 * t^2 + 0.017998 * t^3; 
theta = 2004.3109 * t - 0.42655 * t^2 - 0.041833 * t^3; 
z     = 2306.2181 * t + 1.09468 * t^2 + 0.018203 * t^3; 

% convert arcsec --> deg --> rad 
zeta  = zeta/3600 * pi/180; 
theta = theta/3600 * pi/180; 
z     = z/3600 * pi/180; 

% P row 1 coeffs 
p11 = cos(zeta)*cos(theta)*cos(z) - sin(zeta)*sin(z); 
p12 = -sin(zeta)*cos(theta)*cos(z) - cos(zeta)*sin(z); 
p13 = -sin(theta)*cos(z); 

% P row 2 coeffs 
p21 = cos(zeta)*cos(theta)*sin(z) + sin(zeta)*cos(z); 
p22 = -sin(zeta)*cos(theta)*sin(z) + cos(zeta)*cos(z); 
p23 = -sin(theta)*sin(z); 

% P row 3 coeffs 
p31 = cos(zeta)*cos(theta); 
p32 = -sin(zeta)*sin(theta); 
p33 = cos(theta); 

% P  = precession of ECF wrt ECI 
P = [ p11, p12, p13 ; 
      p21, p22, p23 ; 
      p31, p32, p33 ]; 
  
end