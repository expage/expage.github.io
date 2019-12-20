function out = LSD_fix(lat,lon,age,w,samplingyr,consts)

% Fix w, Rc, and SPhi to be used by P_mu_LSD and LSDspal.
% syntax : LSD_fix(lat,lon,age,w,consts);

% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes, or 0-360 degrees E)
% age = age for production scaling
% w = gravimetric fractional water content - 0.066 is default typically about 14% volumetric per
%     Fred Phillips. -1 gives default value
% samplingyr = year of sampling
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
% Modified by Jakob Heyman (jakob.heyman@gu.se) 2015-2019
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

% Make the time vector
% Age Relative to t0=2010
tv1 = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

Rc = zeros(1,length(tv1));

% Need solar modulation parameter
SPhi = zeros(size(tv1)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

% Fix w
if w < 0;
    out.w = 0.066; % default gravimetric water content for Sato et al. (2008)
else;
    out.w = w;
end;

% interpolate an M for tv1 > 7000...
temp_M = interp1(consts.t_M,consts.M,tv1(77:end));

% catch for negative longitudes before Rc interpolation
if lon < 0; lon = lon + 360;end;

% Make up the Rc vectors.

% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
[loni,lati,tvi] = meshgrid(lon,lat,tv1(1:76));
Rc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,loni,lati,tvi);

% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

Rc(77:end) = temp_M.*(dd(1)*cosd(lat) + ...
    dd(2)*(cosd(lat)).^2 + ...
    dd(3)*(cosd(lat)).^3 + ...
    dd(4)*(cosd(lat)).^4 + ...
    dd(5)*(cosd(lat)).^5 + ...
    dd(6)*(cosd(lat)).^6); 

% constant RcEst for muon production
out.RcEst = (dd(1)*cosd(lat) + ...
    dd(2)*(cosd(lat)).^2 + ...
    dd(3)*(cosd(lat)).^3 + ...
    dd(4)*(cosd(lat)).^4 + ...
    dd(5)*(cosd(lat)).^5 + ...
    dd(6)*(cosd(lat)).^6);

% adjust tv, tv1, Rc, and SPhi to sampling year
if samplingyr <= 2010;
    clipidx = find(tv1 > 2010-samplingyr,1,'first');
    tv = [0 (tv1(clipidx:end)-2010+samplingyr)];
    tv1 = [(tv1-2010+samplingyr) tv1(end)];
    Rc = [Rc Rc(end)];
    SPhi = [SPhi SPhi(end)];
else; % assume 2010 Rc and SPhi value for all years >2010
    tv = [0 (tv1+samplingyr-2010)];
    tv1 = tv;
    Rc = [Rc(1) Rc];
    SPhi = [SPhi(1) SPhi];
end;

% chop off tv at age
clipidx = find(tv < age,1,'last');
out.tv = [tv(1:clipidx) age];

% shorten Rc and SPhi commensurately
out.Rc = interp1(tv1,Rc,out.tv);
out.SPhi = interp1(tv1,SPhi,out.tv);
