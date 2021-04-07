%%  PARTIAL DERIVATIVES 

global muE RE muS AU muM eE wE J2 J3 J4 
global A m p0 r0_drag H  
X  = sym('X', [7 1]); 
dX = sym('dX', [7 1]); 
CD  = X(7); 

% force column vector 
dX(1:3) = X(4:6); 
et = et_t0; 

%% accel due to point mass - DONE

% due to Earth 
r       = sqrt( X(1)^2 + X(2)^2 + X(3)^2 ); 
dX(4:6) = ( - muE / r^3 ) * X(1:3); 

%% accel due to gravity - NOT DONE 

% accel due to J2 (Levinson) 
a0 = -muE;
a2 = -3*J2*RE^2/2; 
a3 = -J3*RE^3/2; 
a4 = -5*J4*RE^4/8; 

x1 = X(1); 
x2 = X(2); 
x3 = X(3); 

dU1 = a0*x1/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 1 ) + a3/r^3*( 35*x3^3/r^3 - 15*x3/r ) + a4/r^4*( 63*x3^4/r^4 - 42*x3^2/r^2 + 3 ) );  
dU2 = a0*x2/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 1 ) + a3/r^3*( 35*x3^3/r^3 - 15*x3/r ) + a4/r^4*( 63*x3^4/r^4 - 42*x3^2/r^2 + 3 ) ); 
dU3 = a0*x3/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 3 ) + a3/r^3*( 35*x3^3/r^3 - 30*x3/r + 3*r/x3 ) + a4/r^4*( 63*x3^4/r^4 - 70*x3^2/r^2 + 15 ) ); 

g_J2J3J4 = fn.g_J2J3J4(X); 

% accel due to J2, J3, and J4  
dX(4:6) = dX(4:6) + [dU1; dU2; dU3]; 

%% accel due to lunisolar perturbation - DONE 

[a_sun, a_moon, X_ESun, X_EMoon] = fn.lunisolar(et, X); 
dX(4:6) = dX(4:6) + a_sun + a_moon; 

%% accel due to SRP 

a_srp = fn.a_SRP(X, X_ESun); 
dX(4:6) = dX(4:6) + a_srp; 

%% accel due to drag - DONE 

pA       = p0 * exp( -(r - r0_drag)/H ); 
Vrel     = [X(4) + wE * X(2); ... 
            X(5) - wE * X(1); ... 
            X(6)]; 
Vnorm    = norm(Vrel); 
a_drag   = - 1/2 * CD * A/m * pA * Vnorm * Vrel;  
dX(4:6) = dX(4:6) + a_drag; 

%% compute partials 

Amat = jacobian( dX, X );       
Amat; 














