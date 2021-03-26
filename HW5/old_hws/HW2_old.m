% ASE 389 Orbit Determination
% HW 2
% Junette Hsin 

clear; 

positionA = [100 100 600 600]; 
positionB = [100 100 700 900]; 

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
dt = 20; 

% set ode45 params 
rel_tol = 3e-14;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-16; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% INTEGRATE! Point mass and J2 
[t_p, x_p] = ode45(@TwoBod_6states, [0:dt:T], [r0; v0], options); 
[t_J2, x_J2] = ode45(@TwoBod_UJ2, [0:dt:T], [r0; v0], options); 

% ------------------------------------------------------------------------

name = 'Problem 1 - 2-Body EOM Orbit'; 
h = figure('name', name, 'position', positionA + [50 0 0 0]); 

    x = x_p; 
    plot3(x(:,1), x(:,2), x(:,3)); hold on; grid on; 
    x = x_J2; 
    plot3(x(:,1), x(:,2), x(:,3)); hold on; grid on; 
    
    x = x_p; 
    plot3(x(1,1), x(1,2), x(1,3), 'mo')
    plot3(x(end,1), x(end,2), x(end,3), 'cx') 
    
    x = x_J2; 
    plot3(x(1,1), x(1,2), x(1,3), 'mo')
    plot3(x(end,1), x(end,2), x(end,3), 'cx') 
    
    legend('point mass', 'J2', 'start', 'end')
    
    xlabel('x (km)')
    ylabel('y (km)') 
    zlabel('z (km)') 
%     legend('orbit', 'start', 'end')
    
    sgtitle(name) 
save_pdf(h, name); 

%% problem 1b 

clear oe oe_check oe_p Tp_p Tp_pJ2 oe_pJ2 

for i = 1:length(x_p) 
    [oe_p(i, :), Tp_p(i,:)] = rv2oe( x_p(i, :) ); 
    oe_check(i, :) = rv2orb_OG( x_p(i, :) ); 
%     Tp_p(i,:) = perigee_pass(oe_p(i,:), x_p(i,:)); 
end 


for i = 1:length(x_J2)
    oe_J2(i,:) = rv2oe( x_J2(i, :) ); 
    oe_check(i,:) = rv2orb_OG( x_J2(i, :) ); 
    
    Tp_J2(i,:) = perigee_pass(oe_J2(i,:), x_J2(i,:)); 
end 

% Time of perigee passage 
% n = mean motion = sqrt(mu/a^3)
% E = acos( r/a * cos(nu) + e )
% M = mean anomaly = E - e*sin(E)
% Tp = t0 - M/n 

% ------------------------------------------------------------------------

labels = {'a', 'e', 'i', '\omega', '\Omega', 'T_p'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad'}; 
name = 'Problem 1b - Orbital Elements'; 
h = figure('name', name, 'position', positionB + [50 0 0 0]); 
for i = 1:5
    subplot(6,1,i) 
        plot(t_J2, oe_J2(:, i)); hold on; grid on; 
        plot(t_p, oe_p(:, i)); 
        title(labels{i}); 
        ylabel(units{i}); 
        increase_ylim; 
        if i == 1
            legend('with J2', 'point mass'); 
        end 
end 
    subplot(6,1,6) 
        plot(t_J2, Tp_J2); hold on; grid on; 
        plot(t_p, Tp_p); 
        title('T_p') 
        ylabel('s') 
sgtitle(name) 
xlabel('Time (sec)') 
save_pdf(h, name); 

%% Problem 1c: energy 

clear vnorm 
clear h 
clear h_mom 
clear dE 
clear dh 
clear dhnorm 

for i = 1:length(x_J2)
    
    U(:,i) = comp_U(x_J2(i, 1:3)); 
    vnorm(:,i) = sqrt( x_J2(i,4)^2 + x_J2(i,5)^2 + x_J2(i,6)^2 ); 
    E(:,i) = vnorm(:,i)^2 / 2 - U(:,i); 
    a = [ x_J2(i,1), x_J2(i,2), x_J2(i,3) ]; 
    b = [ x_J2(i,4), x_J2(i,5), x_J2(i,6) ]; 
    h_mom(:,i) = cross( a , b ); 
    
    dE(:,i) = E(i) - E(1); 
    dh(:,i) = h_mom(:,i) - h_mom(:,1); 
    dhnorm(:,i) = norm( dh(:,i) ); 
    
end 

% ------------------------------------------------------------------------

name = 'Problem 1c - Delta Specific Energy'; 
h = figure('name', name, 'position', positionA + [50 0 0 0]); 
    plot(t_J2, dE); 
    title( 'dE = E(t) - E(t_0)' ) 
    xlabel('Time (s)') 
    ylabel('km^2/s^2') 
save_pdf(h, name); 
    
    
%% Problem 1d: angular momentum 


name = 'Problem 1d - Delta Angular Momentum'; 
h = figure('name', name, 'position', positionA + [50 0 0 0]); 
    plot(t_J2, dh(3, :)); 
    title( 'dh_k = h_k(t) - h_k(t_0)' ) 
    xlabel('Time (s)') 
    ylabel('km^2/s') 
save_pdf(h, name); 
    

%% Problem 2: DRAG 

global CD A m p0 r0_drag H dtheta 

CD = 2; 
A = 3.6e-6; % --> km^2 
m = 1350; 
p0 = 4e-4; % --> -13 --> -4, multiply 1e9 when going 1/m^3 to 1/km^3  
r0_drag = 7298.145; 
H = 200; 
dtheta = 7.29211585530066e-5; 

% set 1 day period (again just in case) 
T = 60 * 60 * 24; 

% set ode45 params 
rel_tol = 3e-14;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-16; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% INTEGRATE! Point mass and J2 and drag 
[t_J2drag, x_J2drag] = ode45(@TwoBod_UJ2_drag, [0:dt:T], [r0; v0], options); 

% INTEGRATE! Point mass and drag 
[t_drag, x_drag] = ode45(@TwoBod_drag, [0:dt:T], [r0; v0], options); 


%% Problem 2a: specific energy 


clear vnorm 
clear h 
clear h_mom 
clear dE 
clear dh 
clear dhnorm 

for i = 1:length(x_J2drag)
    
    U(:,i) = comp_U(x_J2drag(i, 1:3)); 
    vnorm(:,i) = sqrt( x_J2drag(i,4)^2 + x_J2drag(i,5)^2 + x_J2drag(i,6)^2 ); 
    E(:,i) = vnorm(:,i)^2 / 2 - U(:,i); 
    a = [ x_J2drag(i,1), x_J2drag(i,2), x_J2drag(i,3) ]; 
    b = [ x_J2drag(i,4), x_J2drag(i,5), x_J2drag(i,6) ]; 
    h_mom(:,i) = cross( a , b ); 
    
    dE(:,i) = E(i) - E(1); 
    dh(:,i) = h_mom(:,i) - h_mom(:,1); 
    dhnorm(:,i) = norm( dh(:,i) ); 
    
end 

% ------------------------------------------------------------------------

name = 'Problem 2a - Delta Specific Energy'; 
h = figure('name', name, 'position', positionA + [50 0 0 0]); 
    plot(t_J2drag, dE); 
    title( 'dE = E(t) - E(t_0)' ) 
    xlabel('Time (s)') 
    ylabel('km^2/s^2') 
save_pdf(h, name); 
    
%% Problem 2b: orbital elements 


clear oe_drag oe_J2drag Tp_drag 

for i = 1:length(x_drag) 
    [oe_drag(i, :), Tp_drag(i,:)] = rv2oe( x_drag(i, :) ); 
    oe_check(i, :) = rv2orb_OG( x_p(i, :) ); 
%     Tp_p(i,:) = perigee_pass(oe_p(i,:), x_p(i,:)); 
end 

for i = 1:length(x_J2drag) 
    [oe_J2drag(i, :), Tp_J2drag(i,:)] = rv2oe( x_J2drag(i, :) ); 
    oe_check(i, :) = rv2orb_OG( x_p(i, :) ); 
end 

% ------------------------------------------------------------------------

labels = {'a', 'e', 'i', '\omega', '\Omega', 'T_p'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad'}; 
name = 'Problem 2b - Orbital Elements'; 
h = figure('name', name, 'position', positionB + [50 0 0 0]); 
for i = 1:5
    subplot(6,1,i) 
        plot(t_J2drag, oe_J2drag(:,i)); hold on; grid on;
        plot(t_J2, oe_J2(:, i));  
        plot(t_drag, oe_drag(:, i)); 
        plot(t_p, oe_p(:, i)); 
        title(labels{i}); 
        ylabel(units{i}); 
        increase_ylim; 
        if i == 1
            legend('J2 + drag', 'J2', 'drag', 'point mass'); 
        end 
end 
    subplot(6,1,6) 
        plot(t_J2drag, Tp_J2drag); hold on; grid on; 
        plot(t_J2, Tp_J2); 
        plot(t_drag, Tp_drag); 
        plot(t_p, Tp_p); 
        title('T_p') 
        ylabel('s') 
sgtitle(name) 
xlabel('Time (sec)') 
save_pdf(h, name); 



%% ------------------------------------------------------------------------
% differences in orbital elements 

% 6 combinations: 
% (1) p - J2 
% (2) p - drag 
% (3) p - J2drag 
% (4) J2 - drag 
% (5) J2 - J2drag 
% (6) drag - J2drag 

% add Tp, make your life easier 

oe_p = [oe_p, Tp_p]; 
oe_J2 = [oe_J2, Tp_J2]; 
oe_J2drag = [oe_J2drag, Tp_J2drag]; 
oe_drag = [oe_drag, Tp_drag]; 

p_J2 = oe_p - oe_J2; 
p_drag = oe_p - oe_drag; 
p_J2drag = oe_p - oe_J2drag; 
J2_drag = oe_J2 - oe_drag; 
J2_J2drag = oe_J2 - oe_J2drag; 
J2drag_J2 = oe_J2drag - oe_J2; 
drag_J2drag = oe_drag - oe_J2drag; 

labels = {'a', 'e', 'i', '\omega', '\Omega', 'T_p', 'T_p'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad', 's'}; 
name = 'Problem 2b - Orbital Elements Point Mass Diff'; 
h = figure('name', name, 'position', positionB + [50 0 0 0]); 
for i = 1:6
    subplot(6,1,i) 
        plot(t_J2, p_J2(:,i)); hold on; grid on; 
        plot(t_J2, p_drag(:,i));
        plot(t_J2, p_J2drag(:,i)); 
    
        title(labels{i}); 
        ylabel(units{i}); 
        increase_ylim; 
        if i == 1
            legend('p-J2', 'p-drag', 'p-J2drag'); 
        end 
end 
sgtitle(name) 
xlabel('Time (sec)') 
save_pdf(h, name); 

%% 


labels = {'a', 'e', 'i', '\omega', '\Omega', 'T_p', 'T_p'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad', 's'}; 
name = 'Problem 2b - Orbital Elements J2 and Drag Diff'; 
h = figure('name', name, 'position', positionB + [50 0 0 0]); 
for i = 1:6
    subplot(6,1,i) 
%         plot(t_J2, J2_drag(:,i)); hold on; grid on; 
%         plot(t_J2, -J2_J2drag(:,i)); 
        plot(t_J2, J2drag_J2(:,i)); 
%         plot(t_J2, drag_J2drag(:,i)); 
    
        title(labels{i}); 
        ylabel(units{i}); 
        increase_ylim; 
        if i == 1
%             legend('J2-drag', 'J2-J2drag', 'drag-J2drag'); 
            legend('2BJ2drag - 2BJ2') 
        end 
end 
sgtitle(name) 
xlabel('Time (sec)') 
save_pdf(h, name); 



%%     
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

function Tp = perigee_pass(oe, x) 

global mu 

    % perigee passing 
    a = oe(1); 
    e = oe(2); 
    nu = oe(6); 
%     r = norm([ x_pJ2(i,1) x_pJ2(i,2) x_pJ2(i,3) ]); 
    r = sqrt( x(1)^2 + x(2)^2 + x(3)^2 ); 
        
    n = sqrt(mu/a^3); 
    E = acos( r/a * cos(nu) + e );
    M = E - e*sin(E); 
    Tp = M/n; 

end 


