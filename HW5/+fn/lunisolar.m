function [a_sun, a_moon] = lunisolar(et, X)
% From Born 
% Perturbed Motion Equation 2.3.39 

global muS muM 

% Moon wrt Earth 
target      = 'Moon';
frame       = 'J2000';
observer    = 'Earth';
abcorr      = 'NONE';

% get states --> Earth to Moon 
X_Emoon = fn.spice_state(et, target, frame, abcorr, observer); 
X_Emoon = X_Emoon'; 

% Sun wrt Earth 
target  = 'Sun'; 

% get states --> Earth to Sun 
X_Esun = fn.spice_state(et, target, frame, abcorr, observer); 
X_Esun = X_Esun'; 

% pos vector --> sat to sun 
X_sunsat = X(1:6) - X_Esun;
X_satsun = -X_sunsat; 

% Relative position vector of satellite wrt moon 
X_moonsat = X(1:6) - X_Emoon; 
X_satmoon = -X_moonsat; 

% lunisolar acceleration 
a_sun  = muS * ( X_satsun(1:3)/norm(X_satsun(1:3))^3 - X_Esun(1:3)/norm(X_Esun(1:3))^3 );  
a_moon = muM * ( X_satmoon(1:3)/norm(X_satmoon(1:3))^3 - X_Emoon(1:3)/norm(X_Emoon(1:3))^3 ); 

end



