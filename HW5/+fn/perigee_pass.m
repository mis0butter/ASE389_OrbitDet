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