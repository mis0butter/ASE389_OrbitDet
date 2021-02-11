% ASE 389 Orbit Determination
% HW 1 
% Junette Hsin 


%% Problem 1 

t   = linspace( 0, 20, 101 ); 
t   = 0 : 0.2 : 400; 
A   = 1.34; 
phi = pi/3; 
km  = 1; 

x   = A * cos( sqrt( km )*t + phi ); 

name = 'Problem 1: harmonic oscillator x(t)'; 
figure('name', name); 
    plot(t, x); 
    ylabel('x(t)'); 
    xlabel('t'); 
    title(name); 
    set(gca, 'fontsize', 12); 
    
    
%% Problem 2 

y0 = zeros(2,1); 
y0(1) = A * cos(phi); 
y0(2) = - A * sqrt( km ) * sin(phi); 

reltol = 1e-12; 
abstol = 1e-20; 
myoptions   = odeset( 'RelTol', reltol, 'AbsTol', abstol); 
[t, y]      = ode45( @harmoscillator, t, y0, myoptions, km); 

% numerical - analytical 
error = y(:,1)' - x; 

name = 'Problem 2: x(t) error'; 
figure('name', name); 
    plot(t, error); 
    xlabel('t')
    ylabel('numerical - analytical')
    title({name ;
    sprintf('RelTol = %.2e, AbsTol = %.2e ', reltol, abstol) }); 
    set(gca, 'fontsize', 12); 


%% subfunctions 

function dx = harmoscillator( t, x, km )

    dx    = zeros(2, 1); 
    dx(1) = x(2); 
    dx(2) = -km * x(1); 

end 