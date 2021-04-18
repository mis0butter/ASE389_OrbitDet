function y = G_fn(X, XS)
% X  = satellite coords (SAME COORD)
% XS = station coords (SAME COORD)
 
r_site = [X(1)-XS(1); X(2)-XS(2); X(3)-XS(3)]; 
v_site = [X(4)-XS(4); X(5)-XS(5); X(6)-XS(6)]; 
d      = norm(r_site); 
v      = dot(v_site, r_site/norm(r_site)); 

y      = [d; v]; 

end 
    
    