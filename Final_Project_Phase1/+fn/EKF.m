function [t_XSTM, XSTM, Xstar, Y_prefit, Y_postfit, Pi, Lambda] = ... 
    EKF(XSTM_prev, nX, et_t0, t_prop, options, Amat_fn, Ht_fn, P_prev)

global wE LEO_DATA_Apparent W_KJL W_DGO W_ACB r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

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

% % find index for same time state and observation 
t_Y   = LEO_DATA_Apparent(:,2) + et_t0; % time after initial epoch 
ti_X  = t_XSTM(end); 
i_Y   = find(t_Y == ti_X); 
Yi    = LEO_DATA_Apparent(i_Y, :); 

% get JD time 
JD_UTC = cspice_et2utc(t_XSTM(end), 'J', 10); 
JD_UTC = str2num(extractAfter(JD_UTC, 'JD ')); 

% observation covariance 
if Yi(1) == 1
    R = W_KJL; 
    r_STA_ECEF = r_KJL_ECEF; 
elseif Yi(1) == 2
    R = W_DGO; 
    r_STA_ECEF = r_DGO_ECEF; 
else 
    R = W_ACB; 
    r_STA_ECEF = r_ACB_ECEF; 
end 

% Convert station to ECI frame 
r_STA_ECI  = fn.ECEFtoECI(JD_UTC, r_STA_ECEF); 
v_KJL_ECEF = [0; 0; 0]; 
a_ECEF     = v_KJL_ECEF + cross([ 0 0 wE ]', r_STA_ECEF); 
v_STA_ECI  = fn.ECEFtoECI(JD_UTC, a_ECEF); % Technically wrong. Look in Vallado 
XSi        = [r_STA_ECI; v_STA_ECI]; 

% Time update + process noise 
dt    = t_prop(end) - t_prop(1); 
% Q     = diag( (1e-10)^2 * [1 1 1] ); 
Q     = diag( (10000e-10)^2 * [1 1 1] ); 
Gamma = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
P_noise = Gamma * Q * Gamma'; 
Pi_bar = STMi * P_prev * STMi' + P_noise; 

% Y prefit 
Y_prefit(1:2) = Yi(1:2); 
Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 
    
% Observation-state matrix 
Hti = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

% Gain matrix 
Ki = Pi_bar * Hti' * inv( Hti * Pi_bar * Hti' + R ); 
    
% Obtain y difference 
yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 

% Measurement and reference orbit update 
xhat  = Ki * yi; 
Xstar = Xi + xhat; 
Pi    = ( eye(nX) - Ki * Hti ) * Pi_bar;
% Pi    = ( eye(nX) - Ki * Hti ) * Pi_bar * ( eye(nX) - Ki * Hti )' + Ki * R * Ki'; 

% Innovation (information) covariance 
Lambda = (Hti * Pi * Hti' + R); 

% Y postfit 
Y_postfit(1:2) = Yi(1:2); 
Y_postfit(3:4) = fn.G_fn(Xstar, XSi)'; 

end 







