% HW 4

r_ECEF = [ -28738.3218400000; -30844.0723200000; -6.71800000000000 ];
JD = 2458088.50055556; 

[r_ECI] = fn.ECEFtoECI_r(JD,r_ECEF); 