%----------------------- Begin Code Sequence -----------------------------%
% Purpose:                                                                %
% Convert a given set of Keplerian Orbital Elements to state vectors in   %
% in the ECI frame of reference                                           %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%p                      [1 x N]                         Semilatus Rectum
%                                                       (km)
%
%e                      [1 x N]                         Eccentricity Magnitude
%                                                       (unitless)
%
%i                      [1 x N]                         inclination
%                                                       (radians)
%
%O                      [1 x N]                         Right Ascention of
%                                                       the ascending node
%                                                       (radians)
%
%o                      [1 x N]                         Argument of perigee
%                                                       (radians)
%
%nu                     [1 x N]                         True Anomaly
%                                                       (radians)
%
%truLon                 [1 x N]                         True Longitude
%                                                       (radians)
%
%argLat                 [1 x N]                         Argument of Latitude
%                                                       (radians)
%
%lonPer                 [1 x N]                         Longitude of Periapse
%                                                       (radians)
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
%
% Outputs:
%---------                                                                %
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
% References:
%-------------
%Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
%
% Function Dependencies:
%------------------
%None
%------------------------------------------------------------------       %
% Programed by Darin Koblick  03-04-2012                                  %
%------------------------------------------------------------------       %
function [rv] = orb2rv(oe, oe_extra)

a = oe(1) ; 
e = oe(2) ; 
i = oe(3) ; 
w = oe(4) ; 
O = oe(5) ; 
nu = oe(6) ; 

truLon = oe_extra(1); 
argLat = oe_extra(2); 
lonPer = oe_extra(3); 
p = oe_extra(4); 

if ~exist('mu','var');  mu = 398600.4418; end
%Make all inputs consistent w/ dimensions
p = p(:); e = e(:); i = i(:); O = O(:); w = w(:); nu = nu(:);
if exist('truLon','var')
truLon = truLon(:); argLat = argLat(:); lonPer = lonPer(:);
end
ietol = 1e-8;
idx = e < ietol & mod(i,pi) < ietol;
if any(idx); w(idx) = 0; O(idx) = 0; nu(idx) = truLon(idx); end
idx = e < ietol & mod(i,pi) > ietol;
if any(idx); w(idx) = 0; nu(idx) = argLat(idx); end
idx = e > ietol & mod(i,pi) < ietol;
if any(idx); O(idx) = 0; w(idx) = lonPer(idx); end
%Find rPQW and vPQW
rPQW = cat(2, p.*cos(nu)./(1 +e.*cos(nu)), p.*sin(nu)./(1+e.*cos(nu)), zeros(size(nu)));
vPQW = cat(2, -sqrt(mu./p).*sin(nu), sqrt(mu./p).*(e+cos(nu)), zeros(size(nu)));
%Create Transformation Matrix
PQW2IJK = NaN(3,3,size(p,1));
cO = cos(O); sO = sin(O); cw = cos(w); so = sin(w); ci = cos(i); si = sin(i);
PQW2IJK(1,1,:) = cO.*cw-sO.*so.*ci; 
PQW2IJK(1,2,:) = -cO.*so-sO.*cw.*ci; 
PQW2IJK(1,3,:) = sO.*si;
PQW2IJK(2,1,:) = sO.*cw+cO.*so.*ci; 
PQW2IJK(2,2,:) = -sO.*so+cO.*cw.*ci; 
PQW2IJK(2,3,:) = -cO.*si;
PQW2IJK(3,1,:) = so.*si;            
PQW2IJK(3,2,:) = cw.*si;             
PQW2IJK(3,3,:) = ci;
%Transform rPQW and vPQW to rECI and vECI
r = multiDimMatrixMultiply(PQW2IJK,rPQW)';  
v = multiDimMatrixMultiply(PQW2IJK,vPQW)';
rv = [r, v]; 
end

function c = multiDimMatrixMultiply(a,b)
c = NaN(size(b));
c(:,1) = sum(bsxfun(@times,a(1,:,:),permute(b,[3 2 1])),2);
c(:,2) = sum(bsxfun(@times,a(2,:,:),permute(b,[3 2 1])),2);
c(:,3) = sum(bsxfun(@times,a(3,:,:),permute(b,[3 2 1])),2);
end
