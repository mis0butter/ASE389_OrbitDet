function a_srp = a_SRP(et, X) 
% from Vallado 
% Special Perturbation Techniques Equation 8-45 

global RE Cd Cs A m 

% get states --> Earth to Sun 
target   = 'Sun';
frame    = 'J2000';
observer = 'Earth';
abcorr   = 'NONE';
X_Esun   = fn.spice_state(et, target, frame, abcorr, observer); 
X_Esun   = X_Esun'; 

r1 = X(1:3); 
r2 = X_Esun(1:3); 

% Check if X is numeric or symbolic, then check if LOS is in Earth's shadow  
if isnumeric(X)
    % computing numerical - LOS depends 
    tau_min = (norm(r1)^2 - dot(r1, r2)) / ( norm(r1)^2 + norm(r2)^2 - 2*dot(r1, r2) ); 
    LOS = false; 
    if tau_min < 0 || tau_min > 1
        LOS = true; 
    else
        ctau2 = ( (1-tau_min)*norm(r1)^2 + dot(r1, r2)*tau_min ) / RE^2; 
        if ctau2 >= RE^2 
            LOS = true; 
        end 
    end 
else 
    % computing symbolic - LOS true 
    LOS = true; 
end

% CHECK UNITS. USE KM !!!
p_srp     = 4.57e-6;     % solar pressure per unit area, in N/m^2 
p_srp     = p_srp * 1e6; % 1/m^2 --> 1/km^2 
Cr        = 1;     
theta_inc = 0; 
a_srp     = 0; 

if LOS == true 

    % solar panel SRP 
    n = (X_Esun(1:3) - X(1:3))/norm(X_Esun(1:3) - X(1:3));  % normal to solar panel 
    s = X_Esun(1:3)/norm(X_Esun(1:3));                      % sun vector 
    a_srp = -p_srp * A/m * cos(theta_inc) * (2*( Cd/3 + Cs*cos(theta_inc) )*n + (1-Cs)*s); 

end 

end 