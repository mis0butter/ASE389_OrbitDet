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

% initial state initial guess (M --> KM)
CD  = 1.88; 
   
X0 = [ 6984.45711518852 
       1612.2547582643 
       13.0925904314402 
      -1.67667852227336
       7.26143715396544
       0.259889857225218 
       CD ]; 
   
% X0 = [ 85917.5037930985
%           18943.2889929331
%          -3776.58224351347
%         -0.461676934792153
%         -0.937684484459721
%          -1.13962098759127
%          CD ]; 

% initialize STM 
STM0 = eye(7); 
STM0 = reshape(STM0, [49 1]); 
XSTM0 = [X0; STM0];  

global muE RE muS AU muM eE wE J2 J3 J4 Cd Cs eop_data 
global A m p p0 r0_drag H 

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


%% Derive A matrix 

X  = sym('X', [7 1]); 
dX = fn.EOM(et_t0, X); 

% compute partials 
Amat    = jacobian( dX, X );       
Amat_fn = matlabFunction(Amat); 


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
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    X = XSTM(:, 1:6); 
end 
disp('Pos and Vel end: ')
disp(X(end, 1:6)'); 

XSTM_ref0 = XSTM; 


%% convert station coords ECEF --> ECI frame for all t 

% Station coords. Convert M --> KM 
r_KJL_ECEF = [-6143584  1364250  1033743]' / 1000;  % Kwajalein 
r_DGO_ECEF = [ 1907295  6030810 -817119 ]' / 1000;  % Diego 
r_ACB_ECEF = [ 2390310 -5564341  1994578]' / 1000;  % Arecibo 
v_KJL_ECEF = [0; 0; 0]; 
v_DGO_ECEF = [0; 0; 0]; 
v_ACB_ECEF = [0; 0; 0]; 

r_KJL_ECI = []; 
v_KJL_ECI = []; 
r_DGO_ECI = []; 
v_DGO_ECI = []; 
r_ACB_ECI = []; 
v_ACB_ECI = []; 
% t is really et (ephemeris time) 
for i = 1:length(t)
    
    % get current JD time 
    JD_UTC = cspice_et2utc(t(i), 'J', 10); 
    JD_UTC = str2num(extractAfter(JD_UTC, 'JD ')); 
    
    % Convert to ECI frame, save in array 
    r_KJL_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, r_KJL_ECEF); 
    a_ECEF          = v_KJL_ECEF + cross([ 0 0 wE ]', r_KJL_ECEF); 
    v_KJL_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, a_ECEF); % Technically wrong. Look in Vallado 

    % v_KJL_ECI  = fn.ECEFtoECI(JD, v_KJL_ECEF) + cross([ 0 0 wE]', r_KJL_ECEF); 

    r_DGO_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, r_DGO_ECEF); 
    a_ECEF          = v_KJL_ECEF + cross([ 0 0 wE ]', r_DGO_ECEF); 
    v_DGO_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, a_ECEF); % Technically wrong. Look in Vallado 

    r_ACB_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, r_ACB_ECEF); 
    a_ECEF          = v_KJL_ECEF + cross([ 0 0 wE ]', r_ACB_ECEF); 
    v_ACB_ECI(i,:)  = fn.ECEFtoECI(JD_UTC, a_ECEF); % Technically wrong. Look in Vallado 
    
end 

X_KJL_ECI = [r_KJL_ECI, v_KJL_ECI]; 
X_DGO_ECI = [r_DGO_ECI, v_DGO_ECI]; 
X_ACB_ECI = [r_ACB_ECI, v_ACB_ECI]; 


%% DO EVERYTHING IN ECI FRAME 

X  = sym('X', [7; 1]); 
XS = sym('XS', [7; 1]); 
 
r_site = [X(1)-XS(1); X(2)-XS(2); X(3)-XS(3)]; 
v_site = [X(4)-XS(4); X(5)-XS(5); X(6)-XS(6)]; 
d      = norm(r_site); 
v      = dot(v_site, r_site/norm(r_site)); 

Htmat      = sym(zeros(2,7)); 
Htmat(1,:) = simplify(gradient(d, X)); 
Htmat(2,:) = simplify(gradient(v, X)); 
Ht_fn      = matlabFunction(Htmat); 


%% Batch Least Squares a priori 

% initial covariance 
% 10 km std - position 
% 10 m/s    - velocity 

XSTM = XSTM_ref0; 

% Reset t to start incrementing at 0 
t_XSTM = t - t(1); 

% weighting matrices m --> km, mm --> km 
W_KJL = [(10e-3)^2 0; 0 (0.5e-6)^2]; 
W_DGO = [(5e-3)^2  0; 0 (1e-6)^2]; 
W_ACB = [(10e-3)^2 0; 0 (0.5e-6)^2]; 

disp('Find a priori Lambda0 and N0') 
Lambda0 = zeros(7,7); 
N0      = zeros(7,1); 

% Test - Kwajalein 
[Ycalc_KJL, xhat, Lambda0_KJL, N0] = ... 
    fn.batch_LSQ(Yobs_KJL, t_XSTM, XSTM, Ht_fn, X_KJL_ECI, W_KJL, Lambda0, N0); 
norm(xhat(1:3))

% Test - Diego  
[Ycalc_DGO, xhat, Lambda0_DGO, N0] = ... 
    fn.batch_LSQ(Yobs_DGO, t_XSTM, XSTM, Ht_fn, X_DGO_ECI, W_DGO, Lambda0, N0); 
norm(xhat(1:3))

% Test - Arecibo 
[Ycalc_ACB, xhat, Lambda0_ACB, N0] = ... 
    fn.batch_LSQ(Yobs_ACB, t_XSTM, XSTM, Ht_fn, X_ACB_ECI, W_ACB, Lambda0, N0); 
norm(xhat(1:3))


%% Batch all stations 

% xhat = 100; 
XSTM0(1:7) = X0; 
iter = 0; 
while norm(xhat) > 1
    
    % keep track of iterations 
    iter = iter + 1; 
    sprintf('iter = %d', iter)
    
    % update initial conditions 
    XSTM0(1:7) = XSTM0(1:7) + xhat; 
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    
%     % Test - Kwajalein 
%     [Ycalc_KJL, xhat, Lambda0_KJL, N0] = ... 
%         fn.batch_LSQ(Yobs_KJL, t_XSTM, XSTM, Ht_fn, X_KJL_ECI, W_KJL, Lambda0_KJL, N0); 
% 
%     % Test - Diego  
%     [Ycalc_DGO, xhat, Lambda0_DGO, N0] = ... 
%         fn.batch_LSQ(Yobs_DGO, t_XSTM, XSTM, Ht_fn, X_DGO_ECI, W_DGO, Lambda0_DGO, N0); 

    % Test - Arecibo 
    [Ycalc_ACB, xhat, Lambda0_ACB, N0] = ... 
        fn.batch_LSQ(Yobs_ACB, t_XSTM, XSTM, Ht_fn, X_ACB_ECI, W_ACB, Lambda0_ACB, N0); 
    
    disp('xhat pos') 
    xhat(1:3)
    
    disp('x IC pos') 
    XSTM0(1:3)
    
end 

%%

% Diego Garcia 
xhat = 1; 
XSTM0(1:7) = X0; 
iter = 0; 
while norm(xhat) > 1e-1
    
    % keep track of iterations 
    iter = iter + 1; 
    
    % update initial conditions 
    XSTM0(1:7) = XSTM0(1:7) + xhat; 
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    
    % run batch again 
    [Ycalc_DGO, xhat, Lambda, N] = ... 
        fn.batch_LSQ(1, t_XSTM, XSTM, Ht_fn, r_DGO_ECI, v_DGO_ECI, W_DGO); 
    
end 

% Arecibo
xhat = 1; 
XSTM0(1:7) = X0; 
iter = 0; 
while norm(xhat) > 1e-1
    
    % keep track of iterations 
    iter = iter + 1; 
    
    % update initial conditions 
    XSTM0(1:7) = XSTM0(1:7) + xhat; 
    [t, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn), [epochs(1) : 60 : epochs(end)], XSTM0, options); 
    
    % run batch again 
    [Ycalc_ACB, xhat, Lambda, N] = ... 
        fn.batch_LSQ(1, t_XSTM, XSTM, Ht_fn, r_ACB_ECI, v_ACB_ECI, W_ACB); 
    
end 

% Calculate residuals 
[d_err_KJL, d_rms_KJL, v_err_KJL, v_rms_KJL] = calc_res(Yobs_KJL, Ycalc_KJL); 
% Calculate residuals 
[d_err_DGO, d_rms_DGO, v_err_DGO, v_rms_DGO] = calc_res(Yobs_DGO, Ycalc_DGO); 
% Calculate residuals 
[d_err_ACB, d_rms_ACB, v_err_ACB, v_rms_ACB] = calc_res(Yobs_ACB, Ycalc_ACB); 

ftitle = 'Residuals'; 
figure('name', ftitle); 
    subplot(2,1,1) 
        scatter(t_KJL, d_err_KJL); hold on; grid on; 
        scatter(t_DGO, d_err_DGO); 
        scatter(t_ACB, d_err_ACB); 
        title({'Range Residuals. RMS (km): '; ... 
            sprintf('ATL = %.3g, DGO = %.3g, ACB = %.3g', d_rms_KJL, d_rms_DGO, d_rms_ACB) }); 
        ylabel('km')  
    subplot(2,1,2) 
        scatter(t_KJL, v_err_KJL); hold on; grid on; 
        scatter(t_DGO, v_err_DGO); 
        scatter(t_ACB, v_err_ACB); 
        title({'Range-Rate Residuals. RMS (km/s): '; ... 
            sprintf('ATL = %.3g, DGO = %.3g, ACB = %.3g', v_rms_KJL, v_rms_DGO, v_rms_ACB) }); 
        xlabel('Time (s)') 
        ylabel('km/s') 
        legend('Kwajalein', 'Diego', 'Arecibo', 'color', 'none');

%% Subfunctions 

function [d_err_STA, d_rms_STA, v_err_STA, v_rms_STA] = calc_res(Yobs_STA, Ycalc_STA)

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
v_rms_STA = rms(v_err_STA); 

end 













