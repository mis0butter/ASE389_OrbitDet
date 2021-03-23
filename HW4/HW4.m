% HW 4

eop_data = load('finals_iau1980.txt'); 

r_ECEF  = [ -28738.3218400000; -30844.0723200000; -6.71800000000000 ];
JD      = 2458088.50055556; 

[r_ECI] = fn.ECEFtoECI(eop_data, JD, r_ECEF); 

s_vec = [19165.44514777874 -37549.06140374086 -41.043609948282580]'; 
vdiff = r_ECI - s_vec; 
