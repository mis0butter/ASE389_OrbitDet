% HW 5 
% Junette Hsin 

% initial state 
rv = [  6990077.798814194  ; 
        1617465.311978378  ; 
        7273.72441330686   ; 
       -1675.13972506056   ; 
        22679.810569245355 ; 
        252.688512916741   ];
    
% epoch = 1 Feb 2018, 05:00:00 UTC; 
JD = 2458150.70833; 

% Constants 
muE = 398600.4415;          % Earth Gravitational Parameter (km^3/s^2) 
RE  = 6378.1363;            % Earth Radius (km)
muS = 132712440018;         % Sun’s Gravitational Parameter (km^3/s^2)
AU  = 149597870.7;          % 1 Astronomical Unit (km)
muM = 4902.800066;          % Moon’s Gravitational Parameter (km^3/s^2)
eE  = 0.081819221456;       % Earth’s eccentricity 
wE  = 7.292115146706979e-5; % Earth’s rotational velocity (rad/s)

% Atmospheric drag 
r       = norm(rv(1:3)); 
H       = 88667.0;      % m
rho0    = 3.614e-13;    % kg/m3
r0      = (700 + RE);   % m --> km 
rho     = rho0*exp( -(r-r0)/H ); 

% 