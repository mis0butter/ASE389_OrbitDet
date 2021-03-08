function drv_stm = CRTBP_STM(t, rv_stm)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx] time vector (orbit is Keplerian, doesn't matter) 
%   x = [6x] state vector 
% 
% Outputs 
%   dx = [6x] derivative of state vector 
% ------------------------------------------------------------------------

global mu 

% initialize 
drv = zeros(6, 1);   % force column vector 
drv_stm = zeros(42, 1); 

rv = rv_stm(1:6); 
STM = rv_stm(7:42); 
STM = reshape(STM, [6 6]); 

x   = rv(1); 
y   = rv(2); 
z   = rv(3); 
dx  = rv(4); 
dy  = rv(5); 
dz  = rv(6); 

% r_norm  = sqrt(x^2 + y^2 + z^2); 
r1  = sqrt( (x+mu)^2 + y^2 + z^2 ); 
r2  = sqrt( (x-1+mu)^2 + y^2 + z^2 ); 

C1  = (1-mu)/r1^3; 
C2  = mu/r2^3; 

ddx = 2*dy + x - C1*(x+mu) - C2*(x-1+mu); 
ddy = -2*dx + y - C1*y - C2*y; 
ddz = -C1*z - C2*z; 

drv(1:6) = [dx; dy; dz; ddx; ddy; ddz]; 

%% STM stuff 

u1_x = 1 - (1-mu)*(1/(r1^3) - 3*((x+mu)^2)/(r1^5)) - mu*(1/(r2^3) - 3*((x-(1-mu))^2)/(r2^5));
u2_y = 1 - (1-mu)*(1/(r1)^3 - 3*y^2/r1^5) - mu*(1/r2^3-3*y^2/r2^5);
u3_z = (-1)*(1-mu)*(1/(r1)^3 - 3*z^2/r1^5) - mu*(1/r2^3-3*z^2/r2^5);
u1_y = 3*(1-mu)*y*(x+mu)/r1^5 + 3*mu*y*(x-(1-mu))/r2^5;
u1_z = 3*(1-mu)*z*(x+mu)/r1^5 + 3*mu*z*(x-(1-mu))/r2^5;
u2_z = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;

%equality of mixed partials gives (as all the terms are already partials
%of the potential function);
u3_y = u2_z;
u2_x = u1_y;
u3_x = u1_z;

% K 
K = [ 0, 2, 0;
     -2, 0, 0;
      0, 0, 0];

%Then (as mentioned) G is the matrix of partials
G = [u1_x, u1_y, u1_z;
     u2_x, u2_y, u2_z;
     u3_x, u3_y, u3_z];
     
A = [ zeros(3), eye(3); ... 
    G, K ]; 

dSTM = A * STM; 
dSTM = reshape(dSTM, [36 1]); 

drv_stm = [drv; dSTM];  

end