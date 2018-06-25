function LSD_fix_out = LSD_fix(lat,lon,age,w,consts)

% Fix w, Rc, and SPhi to be used by P_mu_LSD and LSDspal.
% syntax : LSD_fix(lat,lon,age,w,consts);

% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes, or 0-360 degrees E)
% age = age for production scaling
% w = gravimetric fractional water content - 0.066 is default typically about 14% volumetric per
%     Fred Phillips. -1 gives default value
% consts = consts_LSD (LSD constants structure)
%
% Input values as scalars
%
% Based on code written by Greg Balco -- Berkeley Geochronology Center
% balcs@bgc.org
% 
% Modified by Brent Goehring and  Nat Lifton -- Purdue University
% nlifton@purdue.edu, bgoehrin@purdue.edu
%
% Modified by Jakob Heyman (jakob.heyman@gu.se) 2015-2018
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

% Make the time vector
% Age Relative to t0=2010
tv1 = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

LSDRc = zeros(1,length(tv1));

% Need solar modulation parameter
this_SPhi = zeros(size(tv1)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

% Fix w
if w < 0;
    LSD_fix_out.w = 0.066; % default gravimetric water content for Sato et al. (2008)
else;
    LSD_fix_out.w = w;
end;

% interpolate an M for tv1 > 7000...
temp_M = interp1(consts.t_M,consts.M,tv1(77:end));

% catch for negative longitudes before Rc interpolation
if lon < 0; lon = lon + 360;end;

% Make up the Rc vectors.

% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
[loni,lati,tvi] = meshgrid(lon,lat,tv1(1:76));
LSDRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,loni,lati,tvi);

% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

LSDRc(77:end) = temp_M.*(dd(1)*cosd(lat) + ...
    dd(2)*(cosd(lat)).^2 + ...
    dd(3)*(cosd(lat)).^3 + ...
    dd(4)*(cosd(lat)).^4 + ...
    dd(5)*(cosd(lat)).^5 + ...
    dd(6)*(cosd(lat)).^6); 

% constant RcEst for muon production
LSDRcEst = (dd(1)*cosd(lat) + ...
    dd(2)*(cosd(lat)).^2 + ...
    dd(3)*(cosd(lat)).^3 + ...
    dd(4)*(cosd(lat)).^4 + ...
    dd(5)*(cosd(lat)).^5 + ...
    dd(6)*(cosd(lat)).^6);

% Next, chop off tv1
clipindex = find(tv1 <= age, 1, 'last' );
tv = tv1(1:clipindex);
if tv(end) < min(age,1E7);
    tv = [tv min(age,1E7)];
end;
% Now shorten the Rc's commensurately
Rc = interp1(tv1,LSDRc,tv);
SPhi = interp1(tv1,this_SPhi,tv);

% output structure
LSD_fix_out.Rc = Rc;
LSD_fix_out.RcEst = LSDRcEst;
LSD_fix_out.SPhi = SPhi;
LSD_fix_out.tv = tv;
