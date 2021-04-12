function [t_STA, d_err_STA, d_rms_STA, v_err_STA, v_rms_STA] = ... 
    Y_residuals(ID_STA, t_XSTM, XSTM, Ht_STA_fn)

global LEO_DATA_Apparent 

% Station data 
i_STA     = find(LEO_DATA_Apparent(:, 1) == ID_STA); 
Yobs_STA  = LEO_DATA_Apparent(i_STA, :); 
t_STA     = Yobs_STA(:, 2); 

% Calculate Y = H * x 
Ycalc_STA = zeros(size(Yobs_STA)); 
Ycalc_STA(:, 1:2) = Yobs_STA(:, 1:2); 
for i = 1:length(i_STA) 
    
    % find t index 
    ti    = Yobs_STA(i,2); 
    i_XSTM = find(t_XSTM == ti); 
    
    % Extract states 
    Xi   = XSTM( i_XSTM, 1:7)'; 
    STMi = XSTM( i_XSTM, 8:7+49 ); 
    STMi = reshape(STMi, [7 7]); 
    
    % compute H 
    Hti = Ht_STA_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6)) * STMi; 
    
    % Y = H * x
    Ycalc_STA(i,3:4) = Hti * Xi; 

end 

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
v_rms_STA = rms(v_err_STA); 

end 