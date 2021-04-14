function Ht_fn_out = Ht_fn(XS)
% XS = station coords 

X = sym('X', [7; 1]); 
 
r_site = [X(1)-XS(1); X(2)-XS(2); X(3)-XS(3)]; 
v_site = [X(4)-XS(4); X(5)-XS(5); X(6)-XS(6)]; 
d      = norm(r_site); 
v      = dot(v_site, r_site/norm(r_site)); 

Htmat      = sym(zeros(2,7)); 
Htmat(1,:) = simplify(gradient(d, X)); 
Htmat(2,:) = simplify(gradient(v, X)); 
Ht_fn_out  = matlabFunction(Htmat); 

end 
    
    