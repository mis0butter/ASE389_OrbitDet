% HW 5 
% Junette Hsin 

clear; clc 
addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

% Load SPICE kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

%% Parameters 

% initial state initial guess (km) 
CD  = 1.88; 
X0 = [ 6990.077798814194  ; 
       1617.465311978378  ; 
       22.679810569245355 ; 
      -1.67513972506056   ; 
       7.27372441330686   ; 
       0.252688512916741  ; 
       CD ];

% initialize STM 
STM0 = eye(7); 
STM0 = reshape(STM0, [49 1]);
XSTM0 = [X0; 0; STM0];  

global muE RE muS AU muM eE wE J2 J3 J4 Cd Cs 
global A m p p0 r0_drag H 
global Cnm Snm 
    
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
m   = 2000;                 % satellite mass (kg) 
Cd = 0.04;                  % diffuse reflection 
Cs = 0.04;                  % specular reflection 

J2 = 1.0826300e-3; 
J3 = -2.5321531e-6; 
J4 = -1.6109876e-6; 

% Atmospheric drag 
r   = norm(X0(1:3)); 
H   = 88667.0;      % m
r0_drag  = (700 + RE);   % m --> km 
p0  = 3.614e-13;    % kg/m3
p   = p0*exp( -(r-r0_drag)/H ); 
A   = 15; % m^2 

% spherical harmonics 
Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end

%% Convert t0 to ET, i.e. seconds past J2000, the base time variable for SPICE. function calls.

%  Epoch for initial conditions 
t0      = 'Feb 1, 2018, 05:00:00 UTC'; 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time. 
et_t0   = cspice_str2et( t0 );

% Moon wrt Earth 
target      = 'Moon';
frame       = 'J2000';
observer    = 'Earth';

% Load observation data 
load('LEO_DATA_Apparent.mat') 

% extract observation epochs 
epochs = LEO_DATA_Apparent(:,2); 
epochs = et_t0 + epochs; 

%% A partials or state 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

global run_state 
run_state = 2; 

if run_state == 0
    A_partials 
elseif run_state == 1
    % fuck the STM for now. integrate 
    [t, X] = ode45(@fn.EOM, [epochs(1) epochs(2)], [X0], options); 
elseif run_state == 2
    % Just in case A_partials run earlier - reverts CD 
    CD  = 1.88; 
    % integrate 
    [t, XSTM] = ode45(@fn.EOM_STM, [epochs(1) epochs(2)], [XSTM0], options); 
end 



