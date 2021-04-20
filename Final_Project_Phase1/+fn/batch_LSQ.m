function [Ycalc_STA, Lambda, N] = ... 
    batch_LSQ(Yobs_STA, t_XSTM, XSTM, et_t0, Ht_fn, Lambda0, N0)

global wE W_KJL W_DGO W_ACB 
global r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

% Initialize calculated Y 
Ycalc_STA = zeros(size(Yobs_STA)); 
Ycalc_STA(:, 1:2) = Yobs_STA(:, 1:2); 

% Set up covariance 
nX     = length(N0); 
Lambda = Lambda0; 
N      = N0; 

% for i = 1:length(Yobs_STA)
for i = 1:28
    
    % find index for same time state and observation 
    Yi   = Yobs_STA(i, :); 
    ti   = Yi(2) + et_t0; 
    i_X  = find(t_XSTM == ti); 
    
    % get JD time 
    JD_UTC = cspice_et2utc(ti, 'J', 10); 
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
    
    % Extract states (all in ECI) 
    Xi   = XSTM( i_X, 1:nX)'; 
    STMi = XSTM( i_X, nX+1 : nX+nX^2 ); 
    STMi = reshape(STMi, [nX nX]); 
    
    % compute H [2x7]
    Hi = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)) * STMi; 
    
    % Accumulate observation 
    Ycalc_STA(i,3:4) = Hi * Xi; 
    
    % Obtain y difference 
    yi = Yobs_STA(i,3:4)' - fn.G_fn(Xi, XSi); 

    % Accumulate covariance 
    Lambda = Lambda + Hi' * inv(R) * Hi; 
    N      = N + Hi' * inv(R) * yi; 

end 

end 