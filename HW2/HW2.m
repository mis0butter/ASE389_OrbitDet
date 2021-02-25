% ASE 389 Orbit Determination
% HW 2
% Junette Hsin 

clear; 

%% Problem 1 

syms x y z 
global mu RE J2 

% constants 
mu = 398600.4;      % G * M1 * M2 
RE = 6378.145;      % Earth radius 
J2 = 0.00108248;    % J2 

% testing 
% syms mu RE J2 

% radius 
r = sqrt(x^2 + y^2 + z^2); 

% U point mass 
Up = mu/r; 

% latitude 
phi = asin(z/r); 

% U J2 
UJ2 = -mu/r * J2 * (RE/r)^2 * ( 3/2 * ( sin(phi) )^2 - 1/2 ); 

% gradient 
d_UJ2 = gradient(UJ2, [x y z]); 
d_UJ2 = simplify(d_UJ2); 

% Initial conditions (km)
r0  = [ -2436.45; -2436.45; 6891.037 ]; 
v0  = [ 5.088611; -5.088611; 0 ]; 
rv0 = [r0; v0]; 

% orbital elements (and sanity check) 
oe0 = rv2oe(rv0);
rv_check2 = oe2rv(oe0); 
[oe_check, oe_extra] = rv2orb_OG(rv0);  
[rv_check] = orb2rv_OG(oe_check, oe_extra); 

% 1 day period 
% a = oe(1); 
% T = abs(2 * pi * sqrt(a^3 / mu));        % period 
T = 60 * 60 * 24; 

% set ode45 params 
rel_tol = 1e-14;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-16; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% INTEGRATE! Point mass and J2 
[t_p, x_p] = ode45(@TwoBod_6states, [0 T], [r0; v0], options); 
[t_pJ2, x_pJ2] = ode45(@TwoBod_UJ2, [0 T], [r0; v0], options); 

% ------------------------------------------------------------------------

name = 'Problem 1: 2-Body EOM Orbit'; 
h = figure('name', name); 

    x = x_p; 
    plot3(x(:,1), x(:,2), x(:,3)); hold on; grid on; 
    x = x_pJ2; 
    plot3(x(:,1), x(:,2), x(:,3)); hold on; grid on; 
    
    x = x_p; 
    plot3(x(1,1), x(1,2), x(1,3), 'mo')
    plot3(x(end,1), x(end,2), x(end,3), 'cx') 
    
    x = x_pJ2; 
    plot3(x(1,1), x(1,2), x(1,3), 'mo')
    plot3(x(end,1), x(end,2), x(end,3), 'cx') 
    
    legend('point mass', 'J2', 'start', 'end')
    
    xlabel('x (km)')
    ylabel('y (km)') 
    zlabel('z (km)') 
%     legend('orbit', 'start', 'end')
    
    sgtitle(name) 

%% problem 1b 

clear oe oe_check oe_p 

for i = 1:length(x_pJ2)
    oe_pJ2(i,:) = rv2oe( x_pJ2(i, :) ); 
    oe_check(i,:) = rv2orb_OG( x_pJ2(i, :) ); 
    
    % perigee passing 
    a = oe_pJ2(1); 
    e = oe_pJ2(2); 
    nu = oe_pJ2(6); 
    r = norm([ x_pJ2(i,1) x_pJ2(i,2) x_pJ2(i,3) ]); 
    
    n = sqrt(mu/a^3); 
    E = acos( r/a * cos(nu) + e );
    M = E - e*sin(E); 
    Tp(i,:) = M/n; 
    
end 

for i = 1:length(x_p) 
    oe_p(i, :) = rv2oe( x_p(i, :) ); 
    oe_check(i, :) = rv2orb_OG( x_p(i, :) ); 
end 

% Time of perigee passage 
% n = mean motion = sqrt(mu/a^3)
% E = acos( r/a * cos(nu) + e )
% M = mean anomaly = E - e*sin(E)
% Tp = t0 - M/n 

% ------------------------------------------------------------------------

labels = {'a', 'e', 'i', '\omega', '\Omega', '\nu'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad'}; 
name = 'Problem 1b: Orbital Elements'; 
h = figure('name', name); 
for i = 1:5
    subplot(5,1,i) 
    plot(t_pJ2, oe_pJ2(:, i)); hold on; grid on; 
    plot(t_p, oe_p(:, i)); 
    title(labels{i}); 
    ylabel(units{i}); 
    increase_ylim; 
    if i == 1
        legend('with J2', 'point mass'); 
    end 
end 
sgtitle(name) 
xlabel('Time (sec)') 

%% Problem 1c: energy 

clear vnorm 
clear h 
clear h_mom 
clear dE 
clear dh 
clear dhnorm 

for i = 1:length(x_pJ2)
    
    U(:,i) = comp_U(x_pJ2(i, 1:3)); 
    vnorm(:,i) = sqrt( x_pJ2(i,4)^2 + x_pJ2(i,5)^2 + x_pJ2(i,6)^2 ); 
    E(:,i) = vnorm(:,i)^2 / 2 - U(:,i); 
    a = [ x_pJ2(i,1), x_pJ2(i,2), x_pJ2(i,3) ]; 
    b = [ x_pJ2(i,4), x_pJ2(i,5), x_pJ2(i,6) ]; 
    h_mom(:,i) = cross( a , b ); 
    
    dE(:,i) = E(i) - E(1); 
    dh(:,i) = h_mom(:,i) - h_mom(:,1); 
    dhnorm(:,i) = norm( dh(:,i) ); 
    
end 

% ------------------------------------------------------------------------

name = 'Problem 1c: Delta Specific Energy'; 
h = figure('name', name); 
    plot(t_pJ2, dE); 
    title( 'dE = E(t) - E(t_0)' ) 
    xlabel('Time (s)') 
    ylabel('km^2/s') 
    
    
%% Problem 1d: angular momentum 


name = 'Problem 1d: Delta Angular Momentum'; 
h = figure('name', name); 
    plot(t_pJ2, dh(3, :)); 
    title( 'dh_k = h_k(t) - h_k(t_0)' ) 
    xlabel('Time (s)') 
    ylabel('km^2/s^2') 
    

%% Problem 2



    
%% subfunctions 

function save_pdf(h, name) 

% save as cropped pdf 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,name,'-dpdf','-r0')
    
end 

function increase_ylim

    ylims  = get(gca, 'ylim');
    yd     = ylims(2) - ylims(1); 
    set(gca, 'ylim', [ylims(1) - 0.2*yd, ylims(2) + 0.2*yd  ]); 

end 

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


