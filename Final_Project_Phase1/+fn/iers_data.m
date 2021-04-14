function [xp_rad, yp_rad, dT] = iers_data(eop_data, JD) 
% From Bulletin A: 
% 
% 1         2       3       4       5       6       7       8       
% year      month   day     MJD     xp      dxp     yp      dyp 
% 
% 9         10      11      12      13      14      15      16 
% UT1-UTC   d       LOD     dLOD    dPsi    ddPsi   deps    ddeps 
%           UT1-UTC         

% From Bulletin B: 
% 17        18          19          20                  21      
% PM x asec PM y asec   UT1-UTC (s) dPsi (milli asec)   deps (milli asec) 

%% get xp, yp 

% Modified JD 
MJD = JD - 2400000.5; 

% find day 
mjd_day = floor(MJD); 

% day fraction 
dfrac = (MJD - mjd_day) / 86400; 

% find MJD row 
i_row = find(eop_data(:,4) == mjd_day, 1, 'first');  

% interpolate 
xp1_asec = eop_data(i_row, 5); 
xp2_asec = eop_data(i_row + 1, 5); 
xp_asec  = (xp2_asec - xp1_asec) * dfrac + xp1_asec; 

yp1_asec = eop_data(i_row, 7); 
yp2_asec = eop_data(i_row + 1, 7); 
yp_asec  = (yp2_asec - yp1_asec) * dfrac + yp1_asec; 

xp_rad   = xp_asec / 3600 * pi/180; 
yp_rad   = yp_asec / 3600 * pi/180; 

%% get dT 

% julian day for Jan 1, 2000 12:00 TT: 2451545
% UT1 = JD + dUT1_UTC 
% dT = UT1 - 2451545 

% interpolate 
dUT1 = eop_data(i_row, 9) / 86400; % seconds --> days  
dUT2 = eop_data(i_row+1, 9) / 86400; % seconds --> days  
dUT = (dUT2 - dUT1) * dfrac + dUT1; 

UT1 = JD + dUT; 
dT = UT1 - 2451545; 

end 