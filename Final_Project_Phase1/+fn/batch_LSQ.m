function [Ycalc_STA, xhat, Lambda, N] = ... 
    batch_LSQ(Yobs_STA, t_XSTM, XSTM, Ht_fn, X_STA_ECI, W_STA, Lambda0, N0)

global wE LEO_DATA_Apparent W_KJL W_DGO W_ACB r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

% Initialize calculated Y 
Ycalc_STA = zeros(size(Yobs_STA)); 
Ycalc_STA(:, 1:2) = Yobs_STA(:, 1:2); 

% Set up covariance 
nX     = length(N0); 
Lambda = Lambda0; 
N      = N0; 

for i = 1:length(Yobs_STA)
% for i = 1:length(10)
    
    % find index for same time state and observation 
    ti  = Yobs_STA(i,2); 
    i_X = find(t_XSTM == ti); 
    
    % Extract states (all in ECI) 
    Xi   = XSTM( i_X, 1:nX)'; 
    STMi = XSTM( i_X, nX+1 : nX+nX^2 ); 
    STMi = reshape(STMi, [nX nX]); 
    XSi  = X_STA_ECI( i_X, : ); 
    
    % compute H [2x7]
    Hi = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)) * STMi; 
    
    % Accumulate observation 
    Ycalc_STA(i,3:4) = Hi * Xi; 
    
    % Obtain y difference 
    yi = Yobs_STA(i,3:4)' - fn.G_fn(Xi, XSi); 

    % Accumulate covariance 
    Lambda = Lambda + Hi' * W_STA * Hi; 
    N      = N + Hi' * W_STA * yi; 

end 

end 