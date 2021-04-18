function [Ycalc_STA, xhat, Lambda, N] = ... 
    EKF(Yobs_STA, t_XSTM, XSTM, Ht_fn, X_STA_ECI, W_STA, Lambda0, N0)

% Station data 
t_STA     = Yobs_STA(:, 2); 

% Calculate Y = H * x 
Ycalc_STA = zeros(size(Yobs_STA)); 
Ycalc_STA(:, 1:2) = Yobs_STA(:, 1:2); 

% Set up lambda 
Lambda = Lambda0; 
N      = N0; 
N      = zeros(7,1); 
    
% find t index 
ti  = Yobs_STA(1,2); 
i_X = find(t_XSTM == ti); 

for i = 1:length(Yobs_STA)
% for i = 1:length(80)
    
    % find t index 
    ti  = Yobs_STA(i,2); 
    i_X = find(t_XSTM == ti); 
    
    % Extract states (all in ECI) 
    Xi   = XSTM( i_X, 1:7)'; 
    STMi = XSTM( i_X, 8:7+49 ); 
    STMi = reshape(STMi, [7 7]); 
    XSi  = X_STA_ECI( i_X, : ); 
    
    % compute H [2x7]
    Hi = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)) * STMi; 
    
    % Accumulate observation 
    Ycalc_STA(i,3:4) = Hi * Xi; 
    
    yi = Yobs_STA(i,3:4)' - fn.G_fn(Xi, XSi); 

    Lambda = Lambda + Hi' * W_STA * Hi; 
    N      = N + Hi' * W_STA * yi; 

end 

% Solve normal equation 
xhat = inv(Lambda) * N; 
% xhat = inv(H_mat'*W_mat*H_mat) * H_mat'*W_mat*y_mat; 

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
v_rms_STA = rms(v_err_STA); 

end 