% ASE 389 Orbit Determination
% HW 2
% Junette Hsin 

%% Problem 1 

syms x y z 
global mu RE J2 

% constants 
mu = 398600.4;      % G * M1 * M2 
RE = 6378.145;      % Earth radius 
J2 = 0.00108248;    % J2 

% testing 
% syms mu RE J2 

% radius 
r = sqrt(x^2 + y^2 + z^2); 

% U point mass 
Up = mu/r; 

% latitude 
phi = asin(z/r); 

% U J2 
UJ2 = -mu/r * J2 * (RE/r)^2 * ( 3/2 * ( sin(phi) )^2 - 1/2 ); 

% gradient 
d_UJ2 = gradient(UJ2, [x y z]); 
d_UJ2 = simplify(d_UJ2); 

% Initial conditions (km)
r  = [ -2436.45; -2436.45; 6891.037 ]; 
v  = [ 5.088611; -5.088611; 0 ]; 
rv = [r; v]; 

% orbital elements (and sanity check) 
% [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb_OG(r, v);  
% [r,v] = orb2rv_OG(p,eMag,i,O,o,nu,truLon,argLat,lonPer,mu); 
oe = rv2oe(rv);
rv_check = oe2rv(oe); 

% 1 day period 
% a = oe(1); 
% T = abs(2 * pi * sqrt(a^3 / mu));        % period 
T = 60 * 60 * 24; 

rel_tol = 1e-14;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-16; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 
[t,x_p] = ode45(@TwoBod_6states, [0 T], [r; v], options); 

[t,x_pJ2] = ode45(@TwoBod_UJ2, [0 T], [r; v], options); 


%% Problem 2





