function Sp = S_rotation(eop_data, JD) 

%% S' = rotation of ECF about angular velocity vector 

[~, ~, dT] = fn.iers_data(eop_data, JD); 
GSMT = 4.894961212823058751375704430 + dT * ... 
    ( 6.300388098984893552276513720 + dT * ... 
    ( 5.075209994113591478053805523e-15 - ... 
    -9.253097568194335640067190688e-24 * dT) ); 
 
aG = GSMT + dpsi * cos(em); 

Sp = [cos(aG), sin(aG), 0; -sin(aG), cos(aG), 0; 0, 0, 1 ]; 

end 