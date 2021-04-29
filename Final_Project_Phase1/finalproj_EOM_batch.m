% HW 5 
% Junette Hsin 

% clear; clc 
addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

% Load SPICE kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

format long g 

%% Parameters 

global muE RE muS AU muM eE wE J2 J3 J4 Cd Cs eop_data 
global A m p p0 r0_drag H CD

% initial state initial guess (M --> KM)
CD  = 1.88;    
% X0 = [ 6984.45711518852 
%        1612.2547582643 
%        13.0925904314402 
%       -1.67667852227336
%        7.26143715396544
%        0.259889857225218 ]; 
   
% ballpark right answer 
% X0 = [  6978.25607108059
%         1616.30079732436
%         19.7187784486054
%        -1.66208620583392
%         7.26104892573667
%         0.270612926043287 ]; 

nX  = length(X0); 

% initialize STM 
STM0 = eye(length(X0)); 
STM0 = reshape(STM0, [length(X0)^2 1]); 
XSTM0 = [X0; STM0];  

% Constants 
muE = 398600.4415;          % Earth Gravitational Parameter (km^3/s^2) 
RE  = 6378.1363;            % Earth Radius (km)
muS = 132712440018;         % Sun’s Gravitational Parameter (km^3/s^2)
AU  = 149597870.7;          % 1 Astronomical Unit (km)
muM = 4902.800066;          % Moon’s Gravitational Parameter (km^3/s^2)
eE  = 0.081819221456;       % Earth’s eccentricity 
wE  = 7.292115146706979e-5; % Earth’s rotational velocity (rad/s)
m   = 2000;                 % satellite mass (kg) 
Cd  = 0.04;                 % diffuse reflection 
Cs  = 0.04;                 % specular reflection 

J2 = 1.08262617385222e-3; 
J3 = -2.53241051856772e-6;
J4 = -1.61989759991697e-6; 

global r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

% Station coords. Convert M --> KM 
r_KJL_ECEF = [-6143584  1364250  1033743]' / 1000;  % Kwajalein 
r_DGO_ECEF = [ 1907295  6030810 -817119 ]' / 1000;  % Diego 
r_ACB_ECEF = [ 2390310 -5564341  1994578]' / 1000;  % Arecibo 

global LEO_DATA_Apparent 
% Load observation data 
load('LEO_DATA_Apparent.mat') 

% Station data 
ID_STA   = 1; 
i_STA    = find(LEO_DATA_Apparent(:, 1) == ID_STA); 
Yobs_KJL = LEO_DATA_Apparent(i_STA, :); 
t_KWJ    = Yobs_KJL(:,2); 

ID_STA   = 2; 
i_STA    = find(LEO_DATA_Apparent(:, 1) == ID_STA); 
Yobs_DGO = LEO_DATA_Apparent(i_STA, :);
t_DGO    = Yobs_DGO(:,2); 

ID_STA   = 3; 
i_STA    = find(LEO_DATA_Apparent(:, 1) == ID_STA); 
Yobs_ACB = LEO_DATA_Apparent(i_STA, :);
t_ACB    = Yobs_DGO(:,2); 

eop_data = load('finals_iau1980.txt'); 

% Atmospheric drag 
r   = norm(X0(1:3));            % km 
H   = 88667.0 / 1000;           % m --> km 
r0_drag  = (700 + RE);          % m --> km 
p0  = 3.614e-13 * 1e9;          % kg/m3 --> kg/km^3 
p   = p0*exp( -(r-r0_drag)/H ); 
A   = 15 / 1e6;                 % km^2 

global Cnm Snm 

% Gravity 
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
% epoch = 23 March 2018, 08:55:03 UTC ; 
t0      = 'March 23, 2018, 08:55:03 UTC'; 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time. 
et_t0   = cspice_str2et( t0 );

% extract observation epochs 
epochs = LEO_DATA_Apparent(:,2); 
epochs = et_t0 + epochs; 


%% Derive A and H matrices 

X  = sym('X', [length(X0) 1]); 
dX = fn.EOM(et_t0, X); 

% compute partials 
Amat    = jacobian( dX, X );       
Amat_fn = matlabFunction(Amat); 

% DO EVERYTHING IN ECI FRAME 

X  = sym('X', [nX; 1]); 
XS = sym('XS', [nX; 1]); 
 
r_site = [X(1)-XS(1); X(2)-XS(2); X(3)-XS(3)]; 
v_site = [X(4)-XS(4); X(5)-XS(5); X(6)-XS(6)]; 
d      = norm(r_site); 
v      = dot(v_site, r_site/norm(r_site)); 

Htmat      = sym(zeros(2,nX)); 
Htmat(1,:) = simplify(gradient(d, X)); 
Htmat(2,:) = simplify(gradient(v, X)); 
Ht_fn      = matlabFunction(Htmat); 

Ht_r_fn  = matlabFunction(Htmat(1,:)); 
Ht_rr_fn = matlabFunction(Htmat(2,:)); 

%% integrate EOM 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
run_state = 2; 
disp('Running sim ...') 

if run_state == 1
    [t, X] = ode45(@fn.EOM, [epochs(1) : 60 : epochs(end)], X0, options); 
elseif run_state == 2
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    X = XSTM(:, 1:6); 
end 
disp('Pos and Vel end: ')
disp(X(end, 1:6)'); 

XSTM_ref0   = XSTM; 
t_XSTM_ref0 = t; 


%% Setting up filters 

% weighting matrices m --> km, mm --> km 
global R_KJL R_DGO R_ACB 
R_KJL = [(10e-3)^2 0; 0 (0.5e-6)^2]; 
R_DGO = [(5e-3)^2  0; 0 (1e-6)^2]; 
R_ACB = [(10e-3)^2 0; 0 (0.5e-6)^2]; 

% initial covariance error 
% 10 km std - position 
% 10 m/s    - velocity 
P0 = [ 10^2*eye(3), zeros(3); 
       zeros(3),  (10e-3)^2*eye(3) ]; 
Lambda0 = inv(P0); 
global Lambda_KJL0 Lambda_DGO0 Lambda_ACB0
Lambda_KJL0 = Lambda0; 
Lambda_DGO0 = Lambda0; 
Lambda_ACB0 = Lambda0; 

global N_KJL0 N_DGO0 N_ACB0
N0    = Lambda0*X0; 
N_KJL0 = N0; 
N_DGO0 = N0; 
N_ACB0 = N0; 


%% Batch all stations to refine IC 

Xstar  = 100*ones(1,3); 
XSTM   = XSTM_ref0; 
t_XSTM = t_XSTM_ref0; 
XSTM0  = XSTM_ref0(1,:)'; 
iter   = 0; 
N_prev = N0; 
Lambda_prev = Lambda0; 

% Batch first 28 measurements 

while norm(Xstar(1:3)) > 0.1
    
    % keep track of iterations 
    iter = iter + 1; 
    sprintf('iter = %d', iter)

    % Test - All stations 
    [Ycalc_all, Lambda, N] = ... 
        fn.batch_LSQ(LEO_DATA_Apparent, t_XSTM, XSTM, et_t0, Ht_fn, Lambda_prev, N_prev); 

%     % Test - Kwajalein 
%     [Ycalc_all, Lambda, N] = ... 
%         fn.batch_LSQ(Yobs_KJL, t_XSTM, XSTM, et_t0, Ht_fn, Lambda0, N0); 
    
    % Solve normal equation 
    Xstar = inv(Lambda) * N; 
    
    % update covariance 
    Lambda_prev = Lambda0; 
    N_prev      = N0; 

    % update initial conditions 
    XSTM0(1:nX) = XSTM0(1:nX) + Xstar; 
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    
    disp('xhat pos'); Xstar(1:3)
    disp('xhat pos norm'); norm(Xstar(1:3))
    disp('x IC pos'); XSTM0(1:3)
    
end 

Lambda_batch = Lambda; 
XSTM_batch   = XSTM; 
XSTM0_batch  = XSTM0; 

