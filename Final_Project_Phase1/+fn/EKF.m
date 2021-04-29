function [t_XSTM, XSTM, X_update, Y_prefit, Y_postfit, P_update, L_pre, L_post] = ... 
    EKF(Yobs_STA, XSTM_prev, nX, et_t0, t_prop, options, Amat_fn, Ht_fn, Ht_r_fn, Ht_rr_fn, P_prev, DATA)

global wE R_KJL R_DGO R_ACB 
global r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

% Integrate ref trajectory and STM from t = i-1 (prev) to t = i (curr) 
if length(t_prop) > 1
    [t_XSTM, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_prop, XSTM_prev, options); 
    Xi   = XSTM(end,1:nX)'; 
    STMi = XSTM(end,nX+1:end); 
else
    t_XSTM = t_prop; 
    XSTM = XSTM_prev; 
    Xi   = XSTM(1:nX); 
    STMi = XSTM(nX+1:end); 
end 
STMi = reshape(STMi, [nX nX]); 

% find index for same time state and observation 
t_Y   = Yobs_STA(:,2) + et_t0; % time after initial epoch 
ti_X  = t_XSTM(end); 
i_Y   = find(t_Y == ti_X); 
Yi    = Yobs_STA(i_Y, :); 

% observation covariance 
if Yi(1) == 1
    R = R_KJL; 
    r_STA_ECEF = r_KJL_ECEF; 
elseif Yi(1) == 2
    R = R_DGO; 
    r_STA_ECEF = r_DGO_ECEF; 
else 
    R = R_ACB; 
    r_STA_ECEF = r_ACB_ECEF; 
end 

% get JD time 
JD_UTC = cspice_et2utc(t_XSTM(end), 'J', 10); 
JD_UTC = str2num(extractAfter(JD_UTC, 'JD ')); 

% Convert station to ECI frame 
XSi = rv_ECEFtoECI(JD_UTC, r_STA_ECEF, [0; 0; 0]); 

% Time update + process noise 
dt    = t_prop(end) - t_prop(1); 
Q     = diag( (1e-10)^2 * [1 1 1] ); 
Gamma = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
P_noise = 0 * Gamma * Q * Gamma'; 
Ppre_bar = STMi * P_prev * STMi' + P_noise; 

% DATA = 0 --> all stations 
% DATA = 1 --> range only 
% DATA = 2 --> range-rate only 
if DATA == 1 
    [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_1(Yi, Xi, XSi, Ht_r_fn, Ppre_bar, R); 
elseif DATA == 2
    [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_2(Yi, Xi, XSi, Ht_rr_fn, Ppre_bar, R); 
else % DATA == 0 
    [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_0(Yi, Xi, XSi, Ht_fn, Ppre_bar, R); 
%     [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_0_ltcorr(Yi, Xi, XSi, Ht_fn, Ppre_bar, R, ... 
%         ti_X, t_prop, XSTM, options, Amat_fn, r_STA_ECEF); 
    
end

end 

%% subfunctions 

function XSi = rv_ECEFtoECI(JD_UTC, r_STA_ECEF, v_STA_ECEF)

global wE 

    r_STA_ECI  = fn.ECEFtoECI(JD_UTC, r_STA_ECEF); 
    a_ECEF     = v_STA_ECEF + cross([ 0 0 wE ]', r_STA_ECEF); 
    v_STA_ECI  = fn.ECEFtoECI(JD_UTC, a_ECEF); % Technically wrong. Look in Vallado p. 228 
    XSi        = [r_STA_ECI; v_STA_ECI]; 

end

function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_0_ltcorr(Yi, Xi, XSi, Ht_fn, Ppre_bar, R, ... 
    ti_X, t_prop, XSTM, options, Amat_fn, r_STA_ECEF)
    
    % Y prefit 
    Y_prefit(1:2) = Yi(1:2); 
    Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 

    % PREFIT Observation-state matrix 
    Hti_pre = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % PREFIT Innovation (information) covariance 
    L_pre = (Hti_pre * Ppre_bar * Hti_pre' + R); 
    
    %%%%%% CORRECT + UPDATE 

    % Light time correction 
    c     = 299792.458; % km/s 
    lt    = Y_prefit(3) / c;  % range / c = delta time (s) 
    t_lt  = ti_X - lt; 
    t_back = [ti_X, t_lt];  
    if length(t_prop) > 1
        XSTM_end = XSTM(end,:); 
    else
        XSTM_end = XSTM; 
    end
    nX = length(Xi); 
    [t_XSTM_corr, XSTM_corr] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_back, XSTM_end, options); 
    Xi_corr   = XSTM_corr(end,1:nX)'; 

    % get JD time 
    JD_UTC = cspice_et2utc(t_XSTM_corr(end), 'J', 10); 
    JD_UTC = str2num(extractAfter(JD_UTC, 'JD ')); 

    % Convert station to ECI frame 
    XSi_corr = rv_ECEFtoECI(JD_UTC, r_STA_ECEF, [0; 0; 0]); 

    % Gain matrix 
%     Ki = Pcorr_bar * Hti_corr' * inv( Hti_corr * Pcorr_bar * Hti_corr' + R ); 
    Ki = Ppre_bar * Hti_pre' * inv( Hti_pre * Ppre_bar * Hti_pre' + R ); 

    % Obtain y difference 
    yi = Yi(3:4)' - fn.G_fn(Xi_corr, XSi); 

    % Measurement and reference orbit update 
    xhat     = Ki * yi; 
    X_update = Xi + xhat; 
    nX       = length(Xi); 
    P_update = ( eye(nX) - Ki * Hti_pre ) * Ppre_bar; 
    % Pi    = ( eye(nX) - Ki * Hti ) * Pi_bar * ( eye(nX) - Ki * Hti )' + Ki * R * Ki'; 

    % POSTFIT Observation-state matrix <-- FIX, USE XSi 
    Hti_post = Ht_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), ... 
        XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % POSTFIT Innovation (information) covariance 
    L_post = (Hti_post * P_update * Hti_post' + R); 

    % Y postfit 
    Y_postfit(1:2) = Yi(1:2); 
    Y_postfit(3:4) = fn.G_fn(X_update, XSi)'; 

end 

function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_0(Yi, Xi, XSi, Ht_fn, Pi_bar, R)

    % Y prefit 
    Y_prefit(1:2) = Yi(1:2); 
    Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 

    % PREFIT Observation-state matrix 
    Hti_pre = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % PREFIT Innovation (information) covariance 
    L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 

    % Gain matrix 
    Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 

    % Obtain y difference 
    yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 

    % Measurement and reference orbit update 
    xhat     = Ki * yi; 
    X_update = Xi + xhat; 
    nX       = length(Xi); 
    P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;

    % POSTFIT Observation-state matrix 
    Hti_post = Ht_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % Innovation (information) covariance 
    L_post = (Hti_post * P_update * Hti_post' + R); 

    % Y postfit 
    Y_postfit(1:2) = Yi(1:2); 
    Y_postfit(3:4) = fn.G_fn(X_update, XSi)'; 

end 

function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_1(Yi, Xi, XSi, Ht_r_fn, Pi_bar, R)

    R = R(1,1); 

    % Y prefit 
    Y_prefit(1:2) = Yi(1:2); 
    y             = fn.G_fn(Xi, XSi)'; 
    Y_prefit(3)   = y(1); 

    % PREFIT Observation-state matrix 
    Hti_pre = Ht_r_fn(Xi(1), Xi(2), Xi(3), XSi(1), XSi(2), XSi(3)); 

    % PREFIT Innovation (information) covariance 
    L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 

    % Gain matrix 
    Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 

    % Obtain y difference (range only) 
    yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 
    yi = yi(1); 

    % Measurement and reference orbit update 
    xhat     = Ki * yi; 
    X_update = Xi + xhat; 
    nX       = length(Xi); 
    P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;

    % POSTFIT Observation-state matrix 
    Hti_post = Ht_r_fn(X_update(1), X_update(2), X_update(3), XSi(1), XSi(2), XSi(3)); 

    % POSTFIT Innovation (information) covariance 
    L_post = (Hti_post * P_update * Hti_post' + R); 

    % Y postfit 
    Y_postfit(1:2) = Yi(1:2); 
    y              = fn.G_fn(X_update, XSi)'; 
    Y_postfit(3)   = y(1); 

end

function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_2(Yi, Xi, XSi, Ht_rr_fn, Pi_bar, R)

    R = R(2,2); 

    % Y prefit 
    Y_prefit(1:2) = Yi(1:2); 
    y             = fn.G_fn(Xi, XSi)'; 
    Y_prefit(3)   = y(2); 

    % PREFIT Observation-state matrix 
    Hti_pre = Ht_rr_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % PREFIT Innovation (information) covariance 
    L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 

    % Gain matrix 
    Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 

    % Obtain y difference (range only) 
    yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 
    yi = yi(2); 

    % Measurement and reference orbit update 
    xhat     = Ki * yi; 
    X_update = Xi + xhat; 
    nX       = length(Xi); 
    P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;

    % POSTFIT Observation-state matrix 
    Hti_post = Ht_rr_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % Innovation (information) covariance 
    L_post = (Hti_post * P_update * Hti_post' + R); 

    % Y postfit 
    Y_postfit(1:2) = Yi(1:2); 
    y              = fn.G_fn(X_update, XSi)'; 
    Y_postfit(3)   = y(2); 

end









