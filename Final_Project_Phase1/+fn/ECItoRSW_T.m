function T = ECItoRSW_T(X_ECI)

r_ECI = X_ECI(1:3); 
v_ECI = X_ECI(4:6); 

R = r_ECI / norm(r_ECI); 
W = cross(r_ECI, v_ECI) / norm(cross(r_ECI, v_ECI)); 
S = cross(W, R); 

T = [R'; S'; W']; 

end 