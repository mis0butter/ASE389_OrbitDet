function U = comp_U(rv) 

global mu J2 RE 

x = rv(1); 
y = rv(2); 
z = rv(3); 

% radius 
r = sqrt(x^2 + y^2 + z^2); 

% U point mass 
Up = mu/r; 

% latitude 
phi = asin(z/r); 

% U J2 
UJ2 = -mu/r * J2 * (RE/r)^2 * ( 3/2 * ( sin(phi) )^2 - 1/2 ); 

% U point mass 
U = Up + UJ2; 

end 