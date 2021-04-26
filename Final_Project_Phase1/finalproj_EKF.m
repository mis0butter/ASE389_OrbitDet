% HW 5 
% Junette Hsin 

% finalproj_EOM_batch; 
close all; 
load X_batch.mat 
load LEO_DATA_Apparent.mat 

addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

% Load SPICE kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

format long g 

% weighting matrices m --> km, mm --> km 
global R_KJL R_DGO R_ACB 
R_KJL = [(10e-3)^2 0; 0 (0.5e-6)^2]; 
R_DGO = [(5e-3)^2  0; 0 (1e-6)^2]; 
R_ACB = [(10e-3)^2 0; 0 (0.5e-6)^2]; 

% SELECT STATIONS AND DATA PARAMS 
STATIONS = 0; 
DATA = 0; 

%% EKF - all observations 

if STATIONS == 0    % use all station data 
    Yobs_STA = LEO_DATA_Apparent; 
elseif STATIONS == 1 
    R_DGO = R_DGO * 1e10; 
    R_ACB = R_ACB * 1e10; 
elseif STATIONS == 2
    R_KJL = R_KJL * 1e10; 
    R_ACB = R_ACB * 1e10; 
else % STATIONS == 3
    R_KJL = R_KJL * 1e10; 
    R_DGO = R_DGO * 1e10; 
end 

et_obs     = Yobs_STA(:,2) + et_t0; 
XSTM_prev  = XSTM0_batch; 
iter       = 0; 
P_prev     = P0; 

X_EKF      = [];     t_X_EKF    = [];     Y_prefit   = [];     Y_postfit  = []; 
Lpre_mat   = [];     Lpost_mat  = [];     sigma3_pre = [];     sigma3_post = []; 

% EKF 
for i = 1:length(et_obs)

    % keep track of iterations 
    iter = iter + 1; 
    sprintf('iter = %d', iter)

    % Propagate state 
    if     i == 1 && et_obs(1) == et_t0; t_prop = et_obs(i); 
    elseif i == 1;                       t_prop = [et_t0 : 60 : et_obs(1) ]; 
    else;                                t_prop = [et_obs(i-1) : 60 : et_obs(i)]; 
    end

    % EKF. All data, range, or range-only 
    [t_XSTM, XSTM, Xstar, Y_pre, Y_post, P, L_pre, L_post] = fn.EKF(Yobs_STA, XSTM_prev, nX, ... 
        epochs(1), t_prop, options, Amat_fn, Ht_fn, Ht_r_fn, Ht_rr_fn, P_prev, DATA); 

    % save states from current iteration 
    if i == 1 && et_obs(1) == et_t0; X_EKF = XSTM(1:nX)'; 
    else;                            X_EKF = [X_EKF; XSTM(:, 1:nX)]; 
    end 

    t_X_EKF   = [t_X_EKF; t_XSTM]; 
    Y_prefit  = [Y_prefit; Y_pre]; 
    Y_postfit = [Y_postfit; Y_post]; 

    % update measurement for next iteration 
    XSTM_prev = [Xstar; STM0]; 
    P_prev    = P;

    % innovations covariance 
    Lpre_mat  = [Lpre_mat; L_pre]; 
    Lpost_mat = [Lpost_mat; L_post]; 

    % 3-sigma STD 
    if DATA == 0
        sigma3_pre  = [sigma3_pre; sqrt(L_pre(1,1))*3, sqrt(L_pre(2,2))*3];  
        sigma3_post = [sigma3_post; sqrt(L_post(1,1))*3, sqrt(L_post(2,2))*3];  
    else  
        sigma3_pre  = [sigma3_pre; sqrt(L_pre(1,1))*3];  
        sigma3_post = [sigma3_post; sqrt(L_post(1,1))*3];  
    end

end 

% Propagate to last period of time 
t_prop = [ et_obs(end) : 60 : et_t0 + 60*60*24 ]; 
[t_XSTM, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_prop, XSTM_prev, options); 
Xi   = XSTM(end,1:nX)'; 
STMi = XSTM(end,nX+1:end); 
STMi = reshape(STMi, [nX nX]); 

% Time update + process noise 
dt       = t_prop(end) - t_prop(1); 
Q        = diag( (1e-10)^2 * [1 1 1] ); 
Gamma    = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
P_noise  = Gamma * Q * Gamma'; 
P_update = STMi * P_prev * STMi' + P_noise; 

% Save propagated states 
X_EKF_prop   = [X_EKF; XSTM(:,1:nX)]; 
t_X_EKF_prop = [t_X_EKF; t_XSTM]; 


%% Plot satellite position 

ftitle = 'JahSat Orbit'; 
figure('name', ftitle); 
    plot3(XSTM_ref0(:,1), XSTM_ref0(:,2), XSTM_ref0(:,3)); hold on; grid on; 
    plot3(XSTM_batch(:,1), XSTM_batch(:,2), XSTM_batch(:,3)); 
    plot3(X_EKF(:,1), X_EKF(:,2), X_EKF(:,3)); hold on; grid on; 
    plot3(X_EKF(1,1), X_EKF(1,2), X_EKF(1,3), 'o')
    xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 
    legend('initial', 'batch', 'EKF', 'prop') 
    title(ftitle)


%% Calculate residuals 

Ypre_KJL = [];      Ypre_DGO = [];      Ypre_ACB = []; 
Ypost_KJL = [];     Ypost_DGO = [];     Ypost_ACB = []; 
sigma3_KJL = [];    sigma3_DGO = [];    sigma3_ACB = []; 

% Extract states that correspond with station measurements 
for i = 1:length(Y_postfit)

    ti  = Y_postfit(i, 2); 
    ti  = ti + et_t0; 
    i_X = find(t_X_EKF == ti); 
    i_STA = Y_postfit(i, 1); 

    if i_STA == 1
        Ypre_KJL   = [Ypre_KJL; Y_prefit(i,:)]; 
        Ypost_KJL  = [Ypost_KJL; Y_postfit(i,:)]; 
        sigma3_KJL = [sigma3_KJL; sigma3_pre(i,:)]; 
    elseif i_STA == 2
        Ypre_DGO   = [Ypre_DGO; Y_prefit(i,:)]; 
        Ypost_DGO  = [Ypost_DGO; Y_postfit(i,:)]; 
        sigma3_DGO = [sigma3_DGO; sigma3_pre(i,:)]; 
    else
        Ypre_ACB   = [Ypre_ACB; Y_prefit(i,:)]; 
        Ypost_ACB  = [Ypost_ACB; Y_postfit(i,:)]; 
        sigma3_ACB = [sigma3_ACB; sigma3_pre(i,:)]; 
    end 

end 

t_KJL = Yobs_KJL(:,2); 
t_DGO = Yobs_DGO(:,2); 
t_ACB = Yobs_ACB(:,2); 
t_ALL = Y_prefit(:,2); 

[dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL] = calc_res_all(Yobs_KJL, Ypre_KJL, DATA); 
[dpre_err_DGO, dpre_rms_DGO, vpre_err_DGO, vpre_rms_DGO] = calc_res_all(Yobs_DGO, Ypre_DGO, DATA); 
[dpre_err_ACB, dpre_rms_ACB, vpre_err_ACB, vpre_rms_ACB] = calc_res_all(Yobs_ACB, Ypre_ACB, DATA); 
[dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL] = calc_res_all(Yobs_KJL, Ypost_KJL, DATA); 
[dpost_err_DGO, dpost_rms_DGO, vpost_err_DGO, vpost_rms_DGO] = calc_res_all(Yobs_DGO, Ypost_DGO, DATA); 
[dpost_err_ACB, dpost_rms_ACB, vpost_err_ACB, vpost_rms_ACB] = calc_res_all(Yobs_ACB, Ypost_ACB, DATA); 
t_KJL = Yobs_KJL(:,2); 
t_DGO = Yobs_DGO(:,2); 
t_ACB = Yobs_ACB(:,2); 

%     if STATIONS == 0

if DATA == 0

    if STATIONS == 0
        ftitle = 'All Stations, All Data Residuals'; 
    elseif STATIONS == 1
        ftitle = 'Kwajalein Only, All Data Residuals'; 
    elseif STATIONS == 2
        ftitle = 'Diego-Garcia Only, All Data Residuals'; 
    else % STATIONS == 3
        ftitle = 'Arecibo Only, All Data Residuals'; 
    end

    % 1 2 3 4 
    % 5 6 7 8 
    % 9 10 11 12 
    % 13 14 15 16 
    figure('name', ftitle); 
        % first row: subplot(4,4,1:4) 
        plot_res_r(4, 4, 1:4, 'PREFIT range residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
            dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, DATA)
        % second row: subplot(4,4,5:8)
        plot_res_rr(4, 4, 5:8, 'PREFIT range-rate residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
            vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, DATA)
        % third row: subplot(4,4,9:12) 
        plot_res_r(4, 4, 9:12, 'POSTFIT range residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, dpost_err_KJL, ... 
            dpost_err_DGO, dpost_err_ACB, dpost_rms_KJL, dpost_rms_DGO, dpost_rms_ACB, DATA)
        % fourth row: subplot(4,4,13:16)
        plot_res_rr(4, 4, 13:16, 'POSTFIT range-rate residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, vpost_err_KJL, ... 
            vpost_err_DGO, vpost_err_ACB, vpost_rms_KJL, vpost_rms_DGO, vpost_rms_ACB, DATA)

        sgtitle(ftitle); 

elseif DATA == 1

    ftitle = 'All Stations, Range Only Residuals'; 
    figure('name', ftitle); 
        % first row: subplot(4,2,1:4) 
        plot_res_r(4, 4, 1:4, 'PREFIT range residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
            dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, DATA)
        % second row: subplot(4,2,5:8) 
        plot_res_r(4, 4, 5:8, 'POSTFIT range residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, dpost_err_KJL, ... 
            dpost_err_DGO, dpost_err_ACB, dpost_rms_KJL, dpost_rms_DGO, dpost_rms_ACB, DATA)

else

    ftitle = 'All Stations, Range-Rate Only Residuals'; 
    figure('name', ftitle); 
        % first row: subplot(4,2,1:4) 
        plot_res_rr(4, 4, 1:4, 'PREFIT range-rate residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
            vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, DATA)
        % second row: subplot(4,2,5:8) 
        plot_res_rr(4, 4, 5:8, 'POSTFIT range-rate residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, vpost_err_KJL, ... 
            vpost_err_DGO, vpost_err_ACB, vpost_rms_KJL, vpost_rms_DGO, vpost_rms_ACB, DATA)

end


%% radial-intrack-cross-track frame transformation for best estimate 

X_best = [-6330.24806665475
           3306.88524247528
           127.818471301279
          -3.43799666089757
          -6.63343176758256
         -0.235596451267133 ]; 

T_best = fn.ECItoRSW_T(X_best); 

dr_ECI = X_best(1:3) - Xi(1:3); 

% transform all measurements 
dr_RSW = T_best * dr_ECI; 

% Plot radial-intrack-crosstrack 
ftitle = 'Radial-Intrack-Crosstrack'; 
figh = figure('name', ftitle); 
    subplot(3,1,1)
        scatter(dr_RSW(1), dr_RSW(2)); hold on; 
        P = P_update(1:2, 1:2); 
        h3 = fn.plot_gaussian_ellipsoid([dr_RSW(1) dr_RSW(2)], P); 
        xlabel('R (km)')
        ylabel('S (km)')
        title('Radial-Intrack') 
    subplot(3,1,2) 
        scatter(dr_RSW(1), dr_RSW(3)); hold on; 
        P = [P_update(1), P_upate(1,3); P_update(3,1), P_update(3,3)]; 
        h3 = fn.plot_gaussian_ellipsoid([dr_RSW(1) dr_RSW(3)], P); 
        xlabel('R (km)')
        ylabel('W (km)')
        title('Radial-Crosstrack') 
    subplot(3,1,3) 
        scatter(dr_RSW(2), dr_RSW(3)); hold on; 
        P = P_update(2:3, 2:3); 
        h3 = fn.plot_gaussian_ellipsoid([dr_RSW(2) dr_RSW(3)], P); 
        xlabel('S (km)')
        ylabel('W (km)')
        title('Intrack-Crosstrack')

        % Vallado ed 4 p. 229 

%% Save data 

if STATIONS == 0
    if DATA == 0
        save('X_ALL_RSW.mat') 
    elseif DATA == 1
        save('X_R_RSW.mat') 
    else
        save('X_RR_RSW.mat') 
    end
elseif STATIONS == 1
    save('X_KJL_RSW.mat') 
elseif STATIONS == 2
    save('X_DGO_RSW.mat') 
else
    save('X_ACB_RSW.mat') 
end 
        
% end

%% Subfunctions 

function plot_res_r(m, n, ivec, ftitle, t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
    dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, DATA)

    subplot(m, n, ivec(1:end-1) )
        scatter(t_KJL, dpre_err_KJL, 3, 'filled'); hold on; grid on; 
        scatter(t_DGO, dpre_err_DGO, 3, 'filled'); 
        scatter(t_ACB, dpre_err_ACB, 3, 'filled'); 
        plot(t_ALL, sigma3_pre(:,1), 'g--'); 
        plot(t_ALL, -sigma3_pre(:,1), 'g--'); 
        title(ftitle) 
        ylabel('km')  
        legend('KJL', 'DGO', 'ACB') 
    subplot(m, n, ivec(end)) 
        stext = { sprintf('KJL mean = %.3g', mean(dpre_err_KJL)); sprintf('KJL RMS = %.3g', dpre_rms_KJL); 
            sprintf('DGO mean = %.3g', mean(dpre_err_DGO)); sprintf('DGO RMS = %.3g', dpre_rms_DGO); 
            sprintf('ACB mean = %.3g', mean(dpre_err_ACB)); sprintf('ACB RMS = %.3g', dpre_rms_ACB)}; 
        text(0, 0.5, stext); axis off

end 

function plot_res_rr(m, n, ivec, ftitle, t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
    vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, DATA)

    if DATA == 0
        sigma3_i = 2; 
    else
        sigma3_i = 1; 
    end

    subplot(m, n, ivec(1:end-1) ) 
        scatter(t_KJL, vpre_err_KJL, 3, 'filled'); hold on; grid on; 
        scatter(t_DGO, vpre_err_DGO, 3, 'filled'); 
        scatter(t_ACB, vpre_err_ACB, 3, 'filled'); 
        plot(t_ALL, sigma3_pre(:,sigma3_i), 'g--'); 
        plot(t_ALL, -sigma3_pre(:,sigma3_i), 'g--'); 
        title(ftitle) 
        ylabel('km')  
        legend('KJL', 'DGO', 'ACB') 
    subplot(m, n, ivec(end) ) 
        stext = { sprintf('KJL mean = %.3g', mean(vpre_err_KJL)); sprintf('KJL RMS = %.3g', vpre_rms_KJL); 
            sprintf('DGO mean = %.3g', mean(vpre_err_DGO)); sprintf('DGO RMS = %.3g', vpre_rms_DGO); 
            sprintf('ACB mean = %.3g', mean(vpre_err_ACB)); sprintf('ACB RMS = %.3g', vpre_rms_ACB)}; 
        text(0, 0.5, stext); axis off

end

function res_all_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL] = calc_res_all(Yobs_KJL, Ypre_KJL); 
    [dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL] = calc_res_all(Yobs_KJL, Ypost_KJL); 

    plot_res_all(ftitle, t_KJL, sigma3_KJL, dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL, ... 
        dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL)

end 

function res_r_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [dpre_err_KJL, dpre_rms_KJL]   = calc_res_r(Yobs_KJL, Ypre_KJL); 
    [dpost_err_KJL, dpost_rms_KJL] = calc_res_r(Yobs_KJL, Ypost_KJL); 

    plot_res_r(ftitle, t_KJL, sigma3_KJL, dpre_err_KJL, dpre_rms_KJL, dpost_err_KJL, dpost_rms_KJL)

end 

function res_rr_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [vpre_err_KJL, vpre_rms_KJL]   = calc_res_rr(Yobs_KJL, Ypre_KJL); 
    [vpost_err_KJL, vpost_rms_KJL] = calc_res_rr(Yobs_KJL, Ypost_KJL); 

    plot_res_r(ftitle, t_KJL, sigma3_KJL, vpre_err_KJL, vpre_rms_KJL, vpost_err_KJL, vpost_rms_KJL)

end 

function [d_err_STA, d_rms_STA, v_err_STA, v_rms_STA] = calc_res_all(Yobs_STA, Ycalc_STA, DATA)

    if DATA == 1 % range only 

        % Calculate residuals 
        d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
        d_rms_STA = rms(d_err_STA); 
        v_err_STA = []; 
        v_rms_STA = []; 

    elseif DATA == 2 % range-rate only 

        % Calculate residuals 
        d_err_STA = []; 
        d_rms_STA = []; 
        v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,3); 
        v_rms_STA = rms(v_err_STA); 

    else

        % Calculate residuals 
        d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
        d_rms_STA = rms(d_err_STA); 
        v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
        v_rms_STA = rms(v_err_STA); 

    end

end 

function [d_err_STA, d_rms_STA] = calc_res_r(Yobs_STA, Ycalc_STA)

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 

end 

function [v_err_STA, v_rms_STA] = calc_res_rr(Yobs_STA, Ycalc_STA)

% Calculate residuals 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,3); 
v_rms_STA = rms(v_err_STA); 

end 

function plot_RSW(dr_RSW, Pi_bar)

end





