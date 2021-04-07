function a_srp = a_SRP(X, X_ESun) 
% from Vallado 
% Special Perturbation Techniques Equation 8-45 

global run_state RE Cd Cs A m 

r1 = X(1:3); 
r2 = X_ESun(1:3); 

if run_state > 0 
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

p_srp     = 4.57e-6;    % solar pressure per unit area, in N/m^2 
Cr        = 1;     
theta_inc = 0; 
a_srp     = 0; 

if LOS == true 

    % solar panel SRP 
    n = X_ESun(1:3)/norm(X_ESun(1:3)); 
    s = X_ESun(1:3)/norm(X_ESun(1:3)); 
    a_srp = -p_srp * A/m * cos(theta_inc) * (2*( Cd/3 + Cs*cos(theta_inc) )*n + (1-Cs)*s); 

end 

end 