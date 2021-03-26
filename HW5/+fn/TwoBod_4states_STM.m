function drvSTM = TwoBod_4states_STM(t, rvSTM)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector (orbit is Keplerian, doesn't matter) 
%   x = [4x1] state vector 
% 
% Outputs 
%   dx = [4x1] derivative of state vector 
% ------------------------------------------------------------------------

global mu 

% initialize 
drv     = zeros(4, 1);       % force column vector 
drvSTM  = zeros(20, 1);   % STM is 4-by-4 --> 16 

STM     = rvSTM(5:20); 
STM     = reshape(STM, [4 4]); 

x  = rvSTM(1); 
y  = rvSTM(2); 
dx = rvSTM(3); 
dy = rvSTM(4); 

% dx1 = x3 
% dx2 = x4 
% dx3 = (-u/r^3) * x1
% dx4 = (-u/r^3) * x2

drv(1:2) = [dx; dy]; 
r        = norm([x, y]); 
drv(3:4) = ( - mu / r^3 ) * [x; y]; 

%% STM stuff 

G = [ -mu/r^3 + 3*mu*x^2/r^5 , 3*mu*x*y/r^5            ; ... 
      3*mu*x*y/r^5           , -mu/r^3 + 3*mu*y^2/r^5  ]; 
  
K = zeros(2,2); 

A = [ zeros(2), eye(2); ... 
      G,        K ]; 
  
dSTM = A * STM; 
dSTM = reshape(dSTM, [16 1]); 

drvSTM = [drv; dSTM]; 

end 